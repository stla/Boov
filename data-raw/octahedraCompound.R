# the compound of five octahedra ####
library(gmp) # we use rational numbers for the vertex coordinates
library(Boov) # to use the 'qsqrtPhi' function
library(rgl)

# all vertices
Vertices <- function(phi, half, O){
  t(cbind(
    c( O , O , -phi )
    , c( O , -phi , O )
    , c( -phi , O , O )
    , c( O , O , phi )
    , c( O , phi , O )
    , c( phi , O , O )
    , c( half , phi*half , phi*phi*half )
    , c( phi*half , phi*phi*half , half )
    , c( phi*phi*half , half , phi*half )
    , c( half , phi*half , -phi*phi*half )
    , c( phi*half , -phi*phi*half , half )
    , c( -phi*phi*half , half , phi*half )
    , c( half , -phi*half , phi*phi*half )
    , c( -phi*half , phi*phi*half , half )
    , c( phi*phi*half , half , -phi*half )
    , c( half , -phi*half , -phi*phi*half )
    , c( -phi*half , -phi*phi*half , half )
    , c( -phi*phi*half , half , -phi*half )
    , c( -half , phi*half , phi*phi*half )
    , c( phi*half , phi*phi*half , -half )
    , c( phi*phi*half , -half , phi*half )
    , c( -half , phi*half , -phi*phi*half )
    , c( phi*half , -phi*phi*half , -half )
    , c( -phi*phi*half , -half , phi*half )
    , c( -half , -phi*half , phi*phi*half )
    , c( -phi*half , phi*phi*half , -half )
    , c( phi*phi*half , -half , -phi*half )
    , c( -half , -phi*half , -phi*phi*half )
    , c( -phi*half , -phi*phi*half , -half )
    , c( -phi*phi*half , -half , -phi*half )
  ))
}
vertices <- Vertices((1+sqrt(5))/2, 0.5, 0)
gmpvertices <- Vertices(qsqrtPhi(17), as.bigq(1L, 2L), as.bigq(0L))
vertices1 <- vertices[1+c(0,1,2,3,4,5), ]
vertices2 <- vertices[1+c(6,10,14,27,25,23), ]
vertices3 <- vertices[1+c(7,11,12,28,26,21), ]
vertices4 <- vertices[1+c(8,9,13,29,24,22), ]
vertices5 <- vertices[1+c(15,16,17,18,19,20), ]
gmpvertices1 <- gmpvertices[1+c(0,1,2,3,4,5), ]
gmpvertices2 <- gmpvertices[1+c(6,10,14,27,25,23), ]
gmpvertices3 <- gmpvertices[1+c(7,11,12,28,26,21), ]
gmpvertices4 <- gmpvertices[1+c(8,9,13,29,24,22), ]
gmpvertices5 <- gmpvertices[1+c(15,16,17,18,19,20), ]

# the face indices are common to all octahedra
faces <- rbind(
  c(1L, 2L, 3L),
  c(1L, 5L, 6L),
  c(6L, 2L, 1L),
  c(1L, 3L, 5L),
  c(4L, 3L, 2L),
  c(6L, 4L, 2L),
  c(3L, 4L, 5L),
  c(6L, 5L, 4L)
)

# define the five octahedra meshes
mesh1 <- list(
  "vertices" = vertices1,
  "gmpvertices" = gmpvertices1,
  "faces" = faces
)
mesh2 <- list(
  "vertices" = vertices2,
  "gmpvertices" = gmpvertices2,
  "faces" = faces
)
mesh3 <- list(
  "vertices" = vertices3,
  "gmpvertices" = gmpvertices3,
  "faces" = faces
)
mesh4 <- list(
  "vertices" = vertices4,
  "gmpvertices" = gmpvertices4,
  "faces" = faces
)
mesh5 <- list(
  "vertices" = vertices5,
  "gmpvertices" = gmpvertices5,
  "faces" = faces
)
meshes <- list(
  mesh1, mesh2, mesh3, mesh4, mesh5
)

octahedraCompound <- list(

  "meshes" = lapply(meshes, function(mesh){
    out <- list("vertices" = mesh[["vertices"]], faces = faces)
    attr(out, "toRGL") <- 3L
    class(out) <- "cgalMesh"
    out
  }),

  "gmpmeshes" = lapply(meshes, function(mesh){
    list("vertices" = mesh[["gmpvertices"]], faces = faces)
  }),

  "rglmeshes" = list(
    tmesh3d(
      "vertices"    = t(vertices1),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices2),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices3),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices4),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices5),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    )
  )
)
