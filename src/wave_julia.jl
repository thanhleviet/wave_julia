module wave_julia

	include("wt.jl");

		include("Transforms.jl");
	#export cwtft
	using Reexport

	@reexport using .WT, .Transforms
end
