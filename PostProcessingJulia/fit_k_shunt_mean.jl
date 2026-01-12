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
lineF_filename = "lineFdata.txt"
lineC_filename = "lineCdata.txt"
lineG_filename = "lineGdata.txt"


###### Relevant routines for line F #############
# parse x_coordinate data for line F
lineF_xcoordinate_data = DataParsing.extract_xcoordinatedata_from_COMSOL(lineF_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical velocity line F
lineF_verticalvelocity_data = DataParsing.extract_verticalvelocity_data_from_COMSOL(lineF_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse vertical acceleration data from line F
lineF_verticalacceleration_data = DataParsing.extract_verticalacceleration_data_from_COMSOL(lineF_filename,t_initial,time_step,t_final,number_of_spatial_points)

# get vertical acceleration at midpoint for F
vF_dot = mean(lineF_verticalacceleration_data,dims=1)
vF_dot = vF_dot'

#### Relevant Routines for line C ############
# parse y_coordinate data for line C
lineC_ycoordinate_data = DataParsing.extract_ycoordinatedata_from_COMSOL(lineC_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse horizontal velocity data from line C (to check if -u_C = u_G)
lineC_horizontalvelocity_data = DataParsing.extract_horizontalvelocity_data_from_COMSOL(lineC_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse pressure from line C
lineC_pressure_data = DataParsing.extract_pressuredata_from_COMSOL(lineC_filename,t_initial,time_step,t_final,number_of_spatial_points)


##### Relevant Routines for line G ##############
# parse y_coordinate data for line G
lineG_ycoordinate_data = DataParsing.extract_ycoordinatedata_from_COMSOL(lineG_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse horizontal velocity data from line G (to check if -u_C = u_G)
lineG_horizontalvelocity_data = DataParsing.extract_horizontalvelocity_data_from_COMSOL(lineG_filename,t_initial,time_step,t_final,number_of_spatial_points)
# parse horizontal acceleration from line G
lineG_pressure_data = DataParsing.extract_pressuredata_from_COMSOL(lineG_filename,t_initial,time_step,t_final,number_of_spatial_points)

########### Plots ###################
# check if -uC = uG
uC = mean(lineC_horizontalvelocity_data,dims=1)
uC = uC'
uG = mean(lineG_horizontalvelocity_data,dims=1)
uG = uG'

uC_and_uG = plot(time_array,-uC,linestyle=:dash,label=:"\$-u_C\$",legend=:outerright)
plot!(time_array,uG,linestyle=:solid,label=:"\$u_G\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$u~(m/sec)\$")
savefig(uC_and_uG,"uC_and_uG_means.pdf")
display(uC_and_uG)

##### Check Equation 22 #######
L = 0.25 # m
Dd = 0.01 # m
shunt_distance = L
rho_water = 1000 # kg/m^3

left_hand_side = vF_dot

pC = mean(lineC_pressure_data,dims=1)
pC = pC'
pG = mean(lineG_pressure_data,dims=1)
pG = pG'
right_hand_side = -(pC.-pG)./(rho_water*shunt_distance)

vFdot_vs_pCpG = plot(time_array,left_hand_side,linestyle=:dash,label=:"\$\\dot{v}_F\$",legend=:outerright)
plot!(time_array,right_hand_side,linestyle=:solid,label=:"\$-\\frac{p_C-p_G}{\\tilde{L}_s\\rho_w}\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$Value~(m/sec^2)\$")
savefig(vFdot_vs_pCpG,"vFdot_vs_pCpG_means.pdf")
display(vFdot_vs_pCpG)

# get sum of squared errors 
vF = mean(lineF_verticalvelocity_data,dims=1)
vF = vF'
difference = left_hand_side .- right_hand_side
fun(k) = sum((difference .+ k[1].*vF).^2)
k0 = rand(1)
result = optimize(fun, k0, NelderMead()) #fminsearch
k1 = Optim.minimizer(result)[1]
k1_vF = k1.*vF

k1_fit_means = plot(time_array,difference,linestyle=:dash,label=:"\$-\\frac{p_C-p_G}{\\tilde{L}_s\\rho_w} + \\dot{v}_F\$",legend=:outerright)
plot!(time_array,-k1_vF,linestyle=:solid,label=:"\$-k_1v_F\$")
xlabel!("\$t~(sec)\$")
ylabel!("\$Value~(m/sec^2)\$")
savefig(k1_fit_means,"k1_fit_means.pdf")
display(k1_fit_means)
