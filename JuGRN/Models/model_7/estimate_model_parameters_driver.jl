# include -
include("Include.jl")

# Setup the timescale of the simulation -
time_start = 0.0
time_stop = 24.0
time_step_size = 0.1
number_of_timesteps = length(time_start:time_step_size:time_stop)

# Calculate Scaled Sensitivity Array
# Load the data dictionary (default parameter values) -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# What is the size of the system?
number_of_states = data_dictionary["number_of_states"]

# main loop -
parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
average_scaled_sensitivity_array = zeros(number_of_states,1)

#=
for (parameter_index,parameter_value) in enumerate(parameter_name_mapping_array)

  # grab the dictionary -
  local_data_dictionary = deepcopy(data_dictionary)

  # Solve the adj simulation -

  (T,X) = AdjSimulation(time_start,time_stop,time_step_size,parameter_index,local_data_dictionary)

  # dump the raw sensitivity arrays to disk -
  # you can modify this to point to some place on disk ...
  data_array = [T X]
  file_path = "./sensitivity/AdjSimulation-P"*string(parameter_index)*".dat"
  writedlm(file_path,data_array)

end
=#

file_path = "./sensitivity"
average_scaled_sensitivity_array = calculate_average_scaled_sensitivity_array(file_path,"AdjSimulation-P",data_dictionary)

# ESTIMATE MODEL PARAMETERS DRIVER
# Setup the simulation timescale -
time_start = 0.0
time_stop = 24.0
time_step_size = 0.1

# Load the data dictionary -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# get my initial parameter guess from the data_dictionary -
initial_parameter_guess = Float64[]
parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
number_of_binding_parameters = length(data_dictionary["binding_parameter_dictionary"])
number_of_control_parameters = length(data_dictionary["control_parameter_dictionary"])

# Update data dictionary to match new parameters before calculating obj
for index = 1:length(parameter_mapping_array)

   parameter_name = parameter_mapping_array[index]

   if index <= number_of_binding_parameters
     push!(initial_parameter_guess,data_dictionary["binding_parameter_dictionary"][parameter_name])
   elseif (index>number_of_binding_parameters && index<=(number_of_binding_parameters+number_of_control_parameters))
     push!(initial_parameter_guess,data_dictionary["control_parameter_dictionary"][parameter_name])
   else
     push!(initial_parameter_guess,data_dictionary[parameter_name])
   end
 end

#=
# get my initial parameter guess from the previous run -
obj_array = readdlm("objective_archive.dat.1")
par_array = readdlm("parameter_archive.dat.1")
min_index = indmin(obj_array)
initial_parameter_guess = par_array[:,min_index]
=#


# Search exposes the *run loop* method -
(objective_archive,parameter_archive) = estimate_model_parameters(objective_function,generation_function,acceptance_function,constraint_function,
  initial_parameter_guess,average_scaled_sensitivity_array; maximum_number_of_iterations=200,show_trace=true)

# convert wrapper to actual array -
number_of_parameters = length(initial_parameter_guess)
parameter_array = zeros(number_of_parameters,1)
for parameter_wrapper in parameter_archive

  parameter_vector = parameter_wrapper.array
  parameter_array = [parameter_array parameter_vector]
end
parameter_array = parameter_array[:,2:end]

# write results to disk -
writedlm("objective_archive.dat.1",objective_archive)
writedlm("parameter_archive.dat.1",parameter_array)


# DETERMINE MODEL ERROR

model_error = objective_function(parameter_array[:,end])
writedlm("model_error.dat",model_error)
