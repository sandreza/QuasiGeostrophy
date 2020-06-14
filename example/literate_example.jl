# # My Example
#
# Here we will try to create a literate example that is also code

# First we import our module
using TrialDocs, Plots

# Now we use the functions in the module
# Specifically the test_args_kw function


# ## Math
# We can test some math
# ```math
#   \frac{1}{2}
# ```
# ## Plotting
# Now we can test plotting
x = randn(3)
scatter(x, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = false, title = "simple graph", xlabel = "horizontal axis", ylabel = "vertical axis")

# !!! tip
#     maybe?
