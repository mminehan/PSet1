function AdjSimulation(time_start,time_stop,time_step_size,parameter_index,data_dictionary)

  # First - run the model to steady-state w/no ATRA forcing -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = [XSS; zeros(size(XSS))]
  data_dictionary["initial_condition_array"] = initial_condition_array;

  (T,X) = SolveAdjBalances(time_start,time_stop,time_step_size,parameter_index,data_dictionary)

  return (T,X)

end
