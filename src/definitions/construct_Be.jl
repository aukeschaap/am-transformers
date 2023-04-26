function construct_Be(c, xs, ys)
    # construct shape function
    Emat = [
        xs(1:3) ys(1:3) [1, 1, 1]
    ] \ UniformScaling(1.);

    # Calculate Bx and By from the solution coefficients and the shape function parameters
    Bx = sum(c .* Emat[2,:]);
    By = -sum(c .* Emat[1,:]);

    return Bx, By
end