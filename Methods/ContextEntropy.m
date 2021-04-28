function  context_entropy = ContextEntropy(p_vector)
context_entropy = 0;
for m = p_vector
    context_entropy = context_entropy - m*log(m);
end
end