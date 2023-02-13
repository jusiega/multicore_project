# i am speed

I have to admit I did not include error handling and I am sorry.

still maybe it was not such a dumb move to devote time to thinking of where errors might be rather than writing hundreds of lines of error handling code which I forgot how to write.

as a result most of the program works (fast) apart from the edge handling inside fblur(). this means for example that unless you like a cropped image the best that you may extract from the bokeh filter is a blur which looks like your camera is bad enough not to catch focus properly and not an artistic kind of effect. naturally for artistic usage it is best to blur just a part of the image which probably is best done in some program with a graphical interface such as gimp. what is interesting about gimp though is that it either doesn't have a bokeh filter or my version of it is too old to include it. either way the bokeh filter is a bit of a failure in this program as it is the only function that really requires the fblur() function and can only be approximated with a separable filter.
