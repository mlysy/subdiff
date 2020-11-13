# generate a fsd model object
model <- fsd_model()
# we can modify the elements of generated object
model$theta_names[1] <- "gamma"
