function r = range_or_one(v)
    r = diff(v);
    if r == 0, r = 1; end
end
