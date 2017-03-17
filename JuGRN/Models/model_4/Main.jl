function estimate_model_parameters(objective_function::Function,generation_function::Function,acceptance_function::Function,constraints_function::Function,
  initial_parameter_guess::Array{Float64,1},average_scaled_sensitivity_array::Array{Float64,2};maximum_number_of_iterations=200,show_trace=true)

  # Initialize -
  iteration_counter = 1                 # setup the calculation -
  objective_archive = Float64[]         # initialize the objective archive -
  parameter_archive = ParameterArray[]  # An array of warppers around the parameter vectors

  # what is my constraint violation? (generate a corrected array)
  initial_parameter_guess = constraints_function(initial_parameter_guess)

  # what is the initial value of the objective function?
  initial_objective_value = objective_function(initial_parameter_guess)
  push!(objective_archive,initial_objective_value)

  # Capture the initial parameters -
  initial_parameter_wrapper = ParameterArray()
  initial_parameter_wrapper.array = initial_parameter_guess
  push!(parameter_archive,initial_parameter_wrapper)

  # current_parameter_guess =
  current_parameter_guess = initial_parameter_guess

  # main loop -
  while (iteration_counter<maximum_number_of_iterations)

    # generate a new parameter guess from the generation_function -
    new_parameter_guess = generation_function(current_parameter_guess,constraints_function)

    # evaluate this performance of this guess -
    objective_function_value = objective_function(new_parameter_guess)

    # do we accept this value of the objective function?
    accept_this_value = acceptance_function(objective_function_value,objective_archive)
    if (accept_this_value == true)

      # We have accepted this value - update the parameter, and objective archive
      parameter_wrapper = ParameterArray()
      parameter_wrapper.array = new_parameter_guess
      push!(parameter_archive,parameter_wrapper)
      push!(objective_archive,objective_function_value)

      # update the current parameter guess w/the best guess so far -
      current_parameter_guess = new_parameter_guess

      # send some stuff to the user -
      if (show_trace == true)
        msg = "Move with performance "*string(objective_function_value)*" has been accepted at iteration "*string(iteration_counter)
        println(msg)
      end
    end

    # update the iteration counter -
    iteration_counter = iteration_counter + 1
  end

  # the parameter, and objective archives get returned to the caller -
  return (objective_archive,parameter_archive)
end
