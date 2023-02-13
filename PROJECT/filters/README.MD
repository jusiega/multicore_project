# i am speed

I have to admit I did not include error handling and I am sorry.

still maybe it was not such a dumb move to devote time to thinking of where errors might be rather than writing hundreds of lines of error handling code which I forgot how to write.

as a result most of the program works (fast) apart from the edge handling inside fblur(). this means for example that unless you like a cropped image the best that you may extract from the bokeh filter is a blur which looks like your camera is bad enough not to catch focus properly and not an artistic kind of effect. naturally for artistic usage it is best to blur just a part of the image which probably is best done in some program with a graphical interface such as gimp. what is interesting about gimp though is that it either doesn't have a bokeh filter or my version of it is too old to include it. either way the bokeh filter is a bit of a failure in this program as it is the only function that really requires the fblur() function and can only be approximated with a separable filter.

the compiled files really should and could have landed here with an earlier timestamp but I can't use git properly. the early timestamp is why I kept the filters.cu file as it is. it is not broken but outputs a misleading line which suggests that there is a default.bmp file in the samples folder, which there isn't now in order to reduce the size of the repository but it can be easily prepared using paint.

this problematic output line can be replaced with "../samples/default.jpg ../output/test.bmp 0 7 0" and the initial value of source variable likewise with "../samples/default.jpg". in practice these initial values might be even impossible to pass on to the rest of the program, I tried it at one point and failed but maybe something else was wrong.

the /a.out is kept as a curiosity and definitely does not work well.

perhaps it would be a good idea to also parallelise the r, g, and b channels filtering. beyond the lazyness that made me not do it there was an idea that this might not be beneficial at all for execution configurations which require any queueing and sequential processing to be done on the gpu so to speak, so for sufficiently large images in essence. perhaps this thought was a late night illusion though.
