# MissionMoon
In this project we aim to optimally control the trajectory of a shuttle to the Moon. Starting at a geostationary orbit around Earth the shuttle should reach the moon while minimizing the used energy and therefore minimizing the additional acceleration that is applied by the shuttle. We describe the formulation of a nonlinear programming (NLP) problem and the reformulations required to deal with the nonlinear equality constraints. For generating solutions to the reformulated NLP the solver CasADi IPOPT (Andersson u. a., 2019) is used.