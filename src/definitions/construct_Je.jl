function construct_je(c, conductivity, source)
    # Calculate eddy current loss
    return  norm(source + conductivity * 1/3 * sum(c));
end