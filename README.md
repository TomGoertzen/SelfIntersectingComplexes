# GAP-Packet Self-Intersections

## Overview of functionalities

- Detection of self-intersections
- Given a complex $X$ with known self-intersections: retriangulate it to obtain a new complex $X'$ that is geometric equivalent without self-intersections
- Compute outer hull and chambers
- Remedying non-manifold parts
- Compute STL files

## Main functions/files

- main.gi: top-level
- detection.gi
- retriangulation.gi
- chambers.gi

- nonmanifold.gi (example_fix_ram_edges)

- mini manual

## Installation

**1.** To get the newest version of this GAP 4 package download the archive file `SelfIntersectingComplexes-x.x.tar.gz` from
>   <https://github.com/TomGoertzen/SelfIntersectingComplexes>

**2.** Locate a `pkg/` directory where GAP searches for packages, see
>   [9.2 GAP Root Directories](https://www.gap-system.org/Manuals/doc/ref/chap9.html#X7A4973627A5DB27D)

in the GAP manual for more information.

**3.** Unpack the archive file in such a `pkg/` directory
which creates a subdirectory called `SelfIntersectingComplexes/`.

**4.** Now you can use the package within GAP by entering `LoadPackage("SelfIntersectingComplexes");` on the GAP prompt.



## Bug reports

Please submit bug reports, feature requests and suggestions via our issue tracker at
>  <https://github.com/TomGoertzen/SelfIntersectingComplexes/issues>

## License

SelfIntersectingComplexes is free software you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. For details, see the file LICENSE distributed as part of this package or see the FSF's own site.

## Acknowledgements

The authors of the package thank Friedrich Rober for helping setting up the GAP-package structure.
