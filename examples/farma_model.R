# generate a farma(1,1) model object
model <- farma_model(1,1)
# we can modify the elements of generated object
model$theta_names[1] <- "gamma"