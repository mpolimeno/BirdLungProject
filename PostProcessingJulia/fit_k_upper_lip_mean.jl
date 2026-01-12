include("DataParsing.jl")
using .DataParsing

using Plots
using Statistics
using Optim
using Random
using LinearAlgebra

t_initial = 0.
time_step = 0.02
t_final = 10.
number_of_spatial_points = 100

time_array = collect(t_initial:time_step:t_final)

# piston's parameters
a0 = 0.00725
f = 0.5
omega = 2*pi*f
zp = a0*cos.(omega*time_array)
zp_dot = -a0*omega*sin.(omega*time_array)
zp_dotdot = -a0*omega^2*cos.(omega*time_array)

# define filenames to be parsed
lineA_filename = "lineAdata.txt"
lineL_filename = "lineLdata.txt"
lineM_filename = "lineMdata.txt"


###### Relevant routines for line A #############
# parse x_coordinate data for line A
lineA_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineA_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical velocity line A
lineA_verticalvelocity_data = DataParsing.extract_verticalvelocity_data_from_COMSOL(lineA_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical acceleration data from line A
lineA_verticalacceleration_data = DataParsing.extract_verticalacceleration_data_from_COMSOL(lineA_filename,t_initial,time_step,t_final,number_of_spatial_points)

# get vertical acceleration at midpoint for F
vA_dot = mean(lineA_verticalacceleration_data,dims=1)
vA_dot = vA_dot'

#### Relevant Routines for line L ############
# parse y_coordinate data for line L
lineL_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineL_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical velocity data from line L (to check if -u_C = u_G)
lineL_verticalvelocity_data = DataParsing.extract_verticalvelocity_data_from_COMSOL(lineL_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse pressure from line L
lineL_pressure_data = DataParsing.extract_pressuredata_from_COMSOL(lineL_filename,t_initial,time_step,t_final,number_of_spatial_points)


##### Relevant Routines for line M ##############
# parse y_coordinate data for line M
lineM_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineM_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical velocity data from line M (to check if -u_C = u_G)
lineM_verticalvelocity_data = DataParsing.extract_verticalvelocity_data_from_COMSOL(lineM_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical acceleration from line M
lineM_pressure_data = DataParsing.extract_pressuredata_from_COMSOL(lineM_filename,t_initial,time_step,t_final,number_of_spatial_points)

########### Plots ###################

##### Check Equation 22 #######
L = 0.25 # m
Dd = 0.01 # m
lineL_ycoordinate_data = DataParsing.extract_ycoordinatedata_from_COMSOL(lineL_filename,t_initial,time_step,t_final,number_of_spatial_points)
lineM_ycoordinate_data = DataParsing.extract_ycoordinatedata_from_COMSOL(lineM_filename,t_initial,time_step,t_final,number_of_spatial_points)
lip_distance = mean(abs.(mean(lineL_ycoordinate_data,dims=1).-mean(lineM_ycoordinate_data,dims=1)))
rho_water = 1000 # kg/m^3

left_hand_side = vA_dot

pL = mean(lineL_pressure_data,dims=1)
pL = pL'
pM = mean(lineM_pressure_data,dims=1)
pM = pM'
right_hand_side = -(pL.-pM)./(rho_water*lip_distance)

vAdot_vs_pLpM = plot(time_array,left_hand_side,linestyle=:dash,label=:"\$\\dot{v}_A\$",legend=:outerright)
plot!(time_array,right_hand_side,linestyle=:solid,label=:"\$-\\frac{p_L-p_M}{\\tilde{\\Delta}_s\\rho_w}\$")
plot!(time_array,-zp_dotdot,linestyle=:dot,label=:"\$\\ddot{z}_p\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$Value~(m/sec^2)\$")
savefig(vAdot_vs_pLpM,"vAdot_vs_pLpM_means.pdf")
display(vAdot_vs_pLpM)

# get sum of squared errors 
vA = mean(lineA_verticalvelocity_data,dims=1)
vA = vA'
difference = left_hand_side .- right_hand_side
fun(k) = sum((difference .+ k[1].*vA).^2)
k0 = rand(1)
result = optimize(fun, k0, NelderMead()) #fminsearch
k1 = Optim.minimizer(result)[1]
k1_vA = k1.*vA

k1_fit_means = plot(time_array,difference,linestyle=:dash,label=:"\$-\\frac{p_L-p_M}{\\tilde{\\Delta}_s\\rho_w} + \\dot{v}_A\$",legend=:outerright)
plot!(time_array,-k1_vA,linestyle=:solid,label=:"\$-k_1v_A\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$Value~(m/sec^2)\$")
savefig(k1_fit_means,"k1_fit_means_LAM.pdf")
display(k1_fit_means)
