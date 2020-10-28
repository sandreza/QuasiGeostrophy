using Plots, GraphRecipes
import Plots: plot
export plot, spectrum

function plot(ϕ::FourierField{S, T}; ylims = ()) where {S, T <: FourierMetaData}
    dims = length(ϕ.metadata.grid.grid)
    if dims == 1
        x = ϕ.metadata.grid.grid[1][:]
        dd = ϕ.metadata.transform.backward * ϕ.data
        plot(x, real.(dd))
        plot!(xlabel = "x")
        plot!(ylabel = ϕ.metadata.name)
        plot!(legend = false, ylims = ylims)
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

plot(ϕ::Field{S, T}; kwargs...) where {S <: FourierField, T} = plot(ϕ.data, kwargs)

function spectrum(ϕ::FourierField)
    dims = length(ϕ.metadata.grid.grid)
    if dims == 1
        f = ϕ.data
        g = log10.(abs.(f))[1:div(length(f),2)+1]
        ymax = maximum(g) * 1.1
        ymin = maximum(g)  - 16
        wvi = collect(1:div(length(f),2)+1) .- 1
        p1 = scatter(wvi, g, label = ϕ.metadata.name)
        scatter!(xlabel = "wavenumber index")
        scatter!(ylabel = "log10(spectral amplitude)")
        scatter!(ylims = (ymin, ymax))
        return p1
    elseif dims == 2
        g = log10.(abs.(ϕ.data))
        cmax = maximum(g) * 1.1
        cmin = maximum(g)  - 16
        clims = (cmin, cmax)
        p1 = heatmap(g, clims = clims, color = :thermometer)
        scatter!(xlabel = "x, wavenumber index")
        scatter!(ylabel = "y, wavenumber index")
        scatter!(title = "log10(spectral amplitude)")
        return p1
    end
    return nothing
end

plot(a::AbstractExpression) = plot(to_expr(a))
plot(a::Equation) = plot(to_expr(a))
