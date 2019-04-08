# generate a fLOC model object
model <- floc_model()
# we can modify the elements of generated object
model$theta_names[1] <- "gamma"