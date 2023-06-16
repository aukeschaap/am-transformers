
using oneAPI
println("using oneAPI.jl version: ", oneAPI.version())

function main()

    for i in 1:2
        f = rand(Float32, 500, 1)
        b = rand(Float32, 500, 1)
        f_gpu = oneArray(f)

        @timev f .+= 1f0
        @timev f_gpu .+= 1f0

        # @timev result = f .* b
        # @timev result_gpu = f_gpu .* b

    end
end


main()
