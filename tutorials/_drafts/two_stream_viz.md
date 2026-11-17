<!-- A draft for a tutorial: the two-stream run in HEXAPIC / VIZ.  It was a
     section of the two-stream example's page; an example's page is now only
     the case's description (the code's examples/).  conf.py excludes
     tutorials/_drafts from the build. -->

## What to look at in HEXAPIC / VIZ

```{figure} ../../_static/images/examples/two_stream_app.png
:name: fig:two-stream-app
:width: 100%
:figclass: hx-shot

Step 300, when the instability saturates. Top left: the wells of the
potential. Top right: the phase space of one beam with the trapped electrons.
Bottom right: the energies, with the first peak of the field energy and the
dip of the kinetic energy.
```

Start with the energies. The field energy grows over the first 300 steps,
reaches a peak and settles. On a logarithmic y axis, which you set with
*Scale* in the view, the rise is close to a straight line. It is less
steep than the fastest mode's growth rate in Results: the field energy also
holds the modes that do not grow. The kinetic energy drops by the amount that the field gains.

Then look at the phase space. Each beam starts as a line at its drift speed.
Where the potential wells form, the line bends into a spiral, because the
electrons that are trapped in a well rotate in it.

The potential shows the wells that trap the electrons. They are uniform
along y.
