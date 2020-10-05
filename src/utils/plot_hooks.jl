using Plots
import Plots: plot
function plot(ϕ::FourierField{S, T}) where {S, T <: FourierMetaData}
    dims = length(ϕ.metadata.grid.grid)
    if dims == 1
        x = ϕ.metadata.grid.grid[1][:]
        dd = ϕ.metadata.transform.backward * ϕ.data
        plot(x, real.(dd))
        plot!(xlabel = "x")
        plot!(ylabel = ϕ.metadata.name)
        plot!(legend = false,)
    elseif dims == 2
        x = ϕ.metadata.grid.grid[1][:]
        y = ϕ.metadata.grid.grid[2][:]
        dd = ϕ.metadata.transform.backward * ϕ.data
        contourf(x, y, real.(dd)', linewidth = 0, color = :thermometer)
        contourf!(xlabel = "x")
        contourf!(ylabel = "y")
        contourf!(title =  ϕ.metadata.name)
    else
        print("Plotting is not supported for fields ")
        print("with dimensions greater ≥ 3")
    end
end

plot(ϕ::Field{S, T}) where {S <: FourierField, T} = plot(ϕ.data)