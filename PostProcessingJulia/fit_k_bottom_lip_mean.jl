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
lineB_filename = "lineBdata.txt"
lineN_filename = "lineNdata.txt"
lineO_filename = "lineOdata.txt"


###### Relevant routines for line B #############
# parse x_coordinate data for line B
lineB_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineB_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical velocity line B
lineB_verticalvelocity_data = DataParsing.extract_verticalvelocity_data_from_COMSOL(lineB_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical acceleration data from line B
lineB_verticalacceleration_data = DataParsing.extract_verticalacceleration_data_from_COMSOL(lineB_filename,t_initial,time_step,t_final,number_of_spatial_points)

# get vertical acceleration at midpoint for B
vB_dot = mean(lineB_verticalacceleration_data,dims=1)
vB_dot = vB_dot'

#### Relevant Routines for line N ############
# parse y_coordinate data for line N
lineN_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineN_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse pressure from line N
lineN_pressure_data = DataParsing.extract_pressuredata_from_COMSOL(lineN_filename,t_initial,time_step,t_final,number_of_spatial_points)


##### Relevant Routines for line M ##############
# parse y_coordinate data for line M
lineO_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineO_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical velocity data from line M (to check if -u_C = u_G)
lineO_verticalvelocity_data = DataParsing.extract_verticalvelocity_data_from_COMSOL(lineO_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical acceleration from line M
lineO_pressure_data = DataParsing.extract_pressuredata_from_COMSOL(lineO_filename,t_initial,time_step,t_final,number_of_spatial_points)

########### Plots ###################

##### Check Equation 22 #######
L = 0.25 # m
Dd = 0.01 # m
lineN_ycoordinate_data = DataParsing.extract_ycoordinatedata_from_COMSOL(lineN_filename,t_initial,time_step,t_final,number_of_spatial_points)
lineO_ycoordinate_data = DataParsing.extract_ycoordinatedata_from_COMSOL(lineO_filename,t_initial,time_step,t_final,number_of_spatial_points)
lip_distance = mean(abs.(mean(lineN_ycoordinate_data,dims=1).-mean(lineO_ycoordinate_data,dims=1)))
rho_water = 1000 # kg/m^3

left_hand_side = vB_dot

pN = mean(lineN_pressure_data,dims=1)
pN = pN'
pO = mean(lineO_pressure_data,dims=1)
pO = pO'
right_hand_side = -(pN.-pO)./(rho_water*lip_distance)

vBdot_vs_pNpO = plot(time_array,left_hand_side,linestyle=:dash,label=:"\$\\dot{v}_B\$",legend=:outerright)
plot!(time_array,right_hand_side,linestyle=:solid,label=:"\$-\\frac{p_N-p_O}{\\tilde{\\Delta}_s\\rho_w}\$")
plot!(time_array,-zp_dotdot,label=:"\$\\ddot{z}_p\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$Value~(m/sec^2)\$")
savefig(vBdot_vs_pNpO,"vBdot_vs_pNpO_means.pdf")
display(vBdot_vs_pNpO)

# get sum of squared errors 
vB = mean(lineB_verticalvelocity_data,dims=1)
vB = vB'
difference = left_hand_side .- right_hand_side
fun(k) = sum((difference .+ k[1].*vB).^2)
k0 = rand(1)
result = optimize(fun, k0, NelderMead()) #fminsearch
k1 = Optim.minimizer(result)[1]
k1_vB = k1.*vB

k1_fit_means = plot(time_array,difference,linestyle=:dash,label=:"\$-\\frac{p_N-p_O}{\\tilde{\\Delta}_s\\rho_w} + \\dot{v}_B\$",legend=:outerright)
plot!(time_array,-k1_vB,linestyle=:solid,label=:"\$-k_1v_B\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$Value~(m/sec^2)\$")
savefig(k1_fit_means,"k1_fit_means_NBO.pdf")
display(k1_fit_means)
