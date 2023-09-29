.. _changelog:


#########
Changelog
#########

******************
0.6.4 (unreleased)
******************

******************
0.6.3 (2023-09-29)
******************

Internal Changes
================
- Add draft GitHub release action. By `Kyle Brindley`_.

******************
0.6.2 (2023-09-20)
******************

Internal Changes
================
-  Clean up conda package CI files after conda-build (:issue:`13`, :merge:`35`). By `Sergio Cordova`_.

******************
0.6.1 (2023-07-24)
******************

Breaking changes
================
- Change project, package, and namespace from 'solver tools' to 'tardigrade solver tools' (:issue:`12`, :merge:`33`). By
  `Kyle Brindley`_.

Internal Changes
================
- Clean up conda-build recipe (:issue:`11`, :merge:`29`). By `Kyle Brindley`_.
- Help CMake find the correct Python executable for conda-build on osx-arm64 (:merge:`30`). By `Kyle Brindley`_.
- Remove compiler as a runtime dependency. The OS-correct standard library package is added as a depedency by
  conda-build (:merge:`31`). By `Kyle Brindley`_.
- Build stdlib variants instead of compiler variants (:merge:`32`). By `Kyle Brindley`_.

******************
0.5.1 (2023-06-20)
******************

Breaking changes
================
- Deploy to the Conda environment preferred ``lib`` directory instead of the CMake linux default ``lib64`` (:issue:`10`,
  :merge:`27`). By `Kyle Brindley`_.

******************
0.4.1 (2023-04-03)
******************

Breaking Changes
================
- Required c++17 (:issues:`8`, :merge:`24`). By `Kyle Brindley`_.

******************
0.3.1 (08-31-2022)
******************

Internal Changes
================
- Build package for multiple compiler versions (:issue: `4`, :merge: `15`). By `Sergio Cordova`_.
- Project configuration and conda build recipe changes to allow macOS builds and conda-build test stage (:merge:`25`).
  By `Kyle Brindley`_.
- Add GCC 11 conda package variant build (:issue:`6`, :merge:`18`). By `Kyle Brindley`_.
- Modifications to the ci environment (:merge:`20`). By `Nathan Miller`_.
- Add GCC 10 conda package variant build (:issue:`7`, :merge:`22`). By `Sergio Cordova`_.

******************
0.3.0 (09-01-2022)
******************

Release
=======
- Released version 0.3.0 (:merge:`12`)

Internal Changes
================
- Build, package, and deploy as a Conda package to the AEA Conda channel (:merge:`9`). By `Nathan Miller`_.
- Added the changelog (:merge:`9`). By `Nathan Miller`_.
- Added the updated environment definition (:merge:`10`). By `Nathan Miller`_.
- Added the updated gitlab-ci.yaml file (:merge:`11`). By `Nathan Miller`_. and `Kyle Brindley`_.
