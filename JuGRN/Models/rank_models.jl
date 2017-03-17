# n = number of models to compare
n = 6

#initialize model_error_array
model_error_array = zeros(n,1)
model_name_array = zeros(n,1)

# build array of model errors
for i = 1:n

  file_path = "./model_"*string(i)*"/model_error.dat"
  model_name_array[i] = i
  model_error_array[i] = readdlm(file_path)[1]

end

ranked_error = sortrows(model_error_array) # assumes Array{Float64,2} input

ranked_models = zeros(n,1)
for j = 1:n
  a= find(model_error_array .== ranked_error[j])
  if length(a) == 1
    ranked_models[j] = a[1]
  elseif length(a) > 1 #handle cases with indistinguishible error
    ranked_models[j:(j+length(a)-1)] = a
  elseif isempty(a) == true
    break
  end
end

print(ranked_models)
