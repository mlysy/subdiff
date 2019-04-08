# generate a fBM model object
model <- fbm_model()
# fit the fBM model
csi_fit(model, dX, dT, Tz, var_calc = TRUE)