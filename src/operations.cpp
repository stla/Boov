#ifndef _BOOVHEADER_
#include "Boov.h"
#endif

template <typename MeshT>
void checkMesh(MeshT mesh, size_t i) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "Mesh n\u00b0" + std::to_string(i) + " self-intersects.";
    Rcpp::stop(msg);
  }
}

template <typename MeshT>
void checkMesh2(MeshT mesh, const std::string& what) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "The " + what + " self-intersects.";
    Rcpp::stop(msg);
  }
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Intersection(const Rcpp::List rmeshes,
                   const bool clean,
                   const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("\u2014 Processing mesh n\u00b01...");
  MeshT mesh_0 = makeSurfMesh<MeshT, PointT>(rmesh, clean, triangulate[0]);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    if(i == 1) {
      checkMesh<MeshT>(meshes[0], 1);
      Message("... done.\n");
    } else {
      checkMesh2<MeshT>(meshes[i - 1], "intersection");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("\u2014 Processing mesh n\u00b0" + meshnum + "...");
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(rmesh_i, clean, triangulate[i]);
    checkMesh<MeshT>(mesh_i, i + 1);
    Message("... done.\n");
    const bool ok = PMP::corefine_and_compute_intersection(
      meshes[i - 1], mesh_i, meshes[i]
    );
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List Intersection_EK(const Rcpp::List rmeshes,
                           const bool clean,
                           const bool normals,
                           const Rcpp::LogicalVector triangulate) {
  EMesh3 mesh =
      Intersection<EK, EMesh3, EPoint3>(rmeshes, clean, triangulate);
  return RSurfTEKMesh(mesh, normals);
}

// [[Rcpp::export]]
Rcpp::List Intersection_Q(const Rcpp::List rmeshes,  // must be triangles
                          const bool clean,
                          const bool normals,
                          const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<QMesh3> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("\u2014 Processing mesh n\u00b01...");
  QMesh3 mesh_0 = makeSurfQMesh(rmesh, clean, triangulate[0]);
  // checkMesh<QMesh3>(mesh_0, 0);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    if(i == 1) {
      checkMesh<QMesh3>(meshes[0], 1);
      Message("... done.\n");
    } else {
      checkMesh2<QMesh3>(meshes[i - 1], "intersection");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("\u2014 Processing mesh n\u00b0" + meshnum + "...");
    QMesh3 mesh_i = makeSurfQMesh(rmesh_i, clean, triangulate[i]);
    checkMesh<QMesh3>(mesh_i, i);
    Message("... done.\n");
    const bool ok = PMP::corefine_and_compute_intersection(
      meshes[i - 1], mesh_i, meshes[i]
    );
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return RSurfTQMesh(meshes[nmeshes - 1], normals);
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Difference(const Rcpp::List rmesh1,
                 const Rcpp::List rmesh2,
                 const bool clean,
                 const bool triangulate1,
                 const bool triangulate2) {
  Message("\u2014 Processing mesh n\u00b01...");
  MeshT smesh1 = makeSurfMesh<MeshT, PointT>(rmesh1, clean, triangulate1);
  checkMesh<MeshT>(smesh1, 1);
  Message("... done.\n");
  Message("\u2014 Processing mesh n\u00b02...");
  MeshT smesh2 = makeSurfMesh<MeshT, PointT>(rmesh2, clean, triangulate2);
  checkMesh<MeshT>(smesh2, 2);
  Message("... done.\n");
  MeshT outmesh;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, outmesh);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return outmesh;
}

// [[Rcpp::export]]
Rcpp::List Difference_EK(const Rcpp::List rmesh1,
                         const Rcpp::List rmesh2,
                         const bool clean,
                         const bool normals,
                         const bool triangulate1,
                         const bool triangulate2) {
  EMesh3 mesh = Difference<EK, EMesh3, EPoint3>(
    rmesh1, rmesh2, clean, triangulate1, triangulate2
  );
  return RSurfTEKMesh(mesh, normals);
}

// [[Rcpp::export]]
Rcpp::List Difference_Q(const Rcpp::List rmesh1,
                        const Rcpp::List rmesh2,
                        const bool clean,
                        const bool normals,
                        const bool triangulate1,
                        const bool triangulate2) {
  Message("\u2014 Processing mesh n\u00b01...");
  QMesh3 smesh1 = makeSurfTQMesh(rmesh1, clean, triangulate1);
  checkMesh<QMesh3>(smesh1, 1);
  Message("... done.\n");
  Message("\u2014 Processing mesh n\u00b02...");
  QMesh3 smesh2 = makeSurfTQMesh(rmesh2, clean, triangulate2);
  checkMesh<QMesh3>(smesh2, 2);
  Message("... done.\n");
  QMesh3 outmesh;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, outmesh);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return RSurfQMesh(outmesh, normals);
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Union(const Rcpp::List rmeshes,
            const bool clean,
            const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("\u2014 Processing mesh n\u00b01...");
  MeshT mesh_0 = makeSurfMesh<MeshT, PointT>(rmesh, clean, triangulate[0]);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    if(i == 1) {
      checkMesh<MeshT>(meshes[0], 1);
      Message("... done.\n");
    } else {
      checkMesh2<MeshT>(meshes[i - 1], "union");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("\u2014 Processing mesh n\u00b0" + meshnum + "...");
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(rmesh_i, clean, triangulate[i]);
    checkMesh<MeshT>(mesh_i, i + 1);
    Message("... done.\n");
    const bool ok =
        PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List Union_EK(const Rcpp::List rmeshes,
                    const bool clean,
                    const bool normals,
                    const Rcpp::LogicalVector triangulate) {
  EMesh3 mesh = Union<EK, EMesh3, EPoint3>(rmeshes, clean, triangulate);
  return RSurfTEKMesh(mesh, normals);
}

// [[Rcpp::export]]
Rcpp::List Union_Q(const Rcpp::List rmeshes,
                   const bool clean,
                   const bool normals,
                   const Rcpp::LogicalVector triangulate) {
  const size_t nmeshes = rmeshes.size();
  std::vector<QMesh3> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Message("\u2014 Processing mesh n\u00b01...");
  QMesh3 mesh_0 = makeSurfQMesh(rmesh, clean, triangulate[0]);
  meshes[0] = mesh_0;
  for(size_t i = 1; i < nmeshes; i++) {
    if(i == 1) {
      checkMesh<QMesh3>(meshes[0], 1);
      Message("... done.\n");
    } else {
      checkMesh2<QMesh3>(meshes[i - 1], "union");
    }
    const std::string meshnum = std::to_string(i + 1);
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Message("\u2014 Processing mesh n\u00b0" + meshnum + "...");
    QMesh3 mesh_i = makeSurfQMesh(rmesh_i, clean, triangulate[i]);
    checkMesh<QMesh3>(mesh_i, i + 1);
    Message("... done.\n");
    const bool ok =
        PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }
  return RSurfTQMesh(meshes[nmeshes - 1], normals);
}
