function testfn(one,two,options)
arguments
    one
    two
    options.three;
    options.four = 4;
end
disp(one*options.four)
%disp(options.three * options.four)
end