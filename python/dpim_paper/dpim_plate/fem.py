from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, cg, splu, spsolve

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def poly_row(x: float, y: float) -> FloatArray:
    return np.array([1.0, x, y, x * x, x * y, y * y, x * x * y, x * y * y])


def dnormal_row(x: float, y: float, nx: float, ny: float) -> FloatArray:
    return np.array(
        [
            0.0,
            nx,
            ny,
            2.0 * x * nx,
            y * nx + x * ny,
            2.0 * y * ny,
            2.0 * x * y * nx + x * x * ny,
            y * y * nx + 2.0 * x * y * ny,
        ]
    )


def bending_b(x: float, y: float, cmat: FloatArray) -> FloatArray:
    rxx = np.array([0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0 * y, 0.0])
    ryy = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0 * x])
    rxy = np.array([0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 4.0 * x, 4.0 * y])
    return np.vstack((rxx @ cmat, ryy @ cmat, rxy @ cmat))


def quad_map(xi: float, eta: float, xe: FloatArray, ye: FloatArray) -> tuple[float, float, float]:
    shape = 0.25 * np.array(
        [
            (1.0 - xi) * (1.0 - eta),
            (1.0 + xi) * (1.0 - eta),
            (1.0 + xi) * (1.0 + eta),
            (1.0 - xi) * (1.0 + eta),
        ]
    )
    dxi = 0.25 * np.array([-(1.0 - eta), 1.0 - eta, 1.0 + eta, -(1.0 + eta)])
    deta = 0.25 * np.array([-(1.0 - xi), -(1.0 + xi), 1.0 + xi, 1.0 - xi])
    xq = float(shape @ xe)
    yq = float(shape @ ye)
    jacobian = np.array(
        [[dxi @ xe, deta @ xe], [dxi @ ye, deta @ ye]], dtype=np.float64
    )
    return xq, yq, float(abs(np.linalg.det(jacobian)))


@dataclass(frozen=True, slots=True)
class FEMSolveAudit:
    mean_iterations: float
    max_iterations: int
    fallback_count: int
    max_relative_residual: float
    direct_validation_max_relative_error: float


@dataclass(slots=True)
class PlateModel:
    nodes: FloatArray
    quads: IntArray
    n_free: int
    free_dofs: IntArray
    force: FloatArray
    response_vector: FloatArray
    csr_indices: IntArray
    csr_indptr: IntArray
    assembly_slot: IntArray
    contribution_element: IntArray
    contribution_base_value: FloatArray
    mean_data: FloatArray
    mean_matrix: csr_matrix
    mean_solution: FloatArray
    _mean_lu: object

    @property
    def nnz(self) -> int:
        return int(self.mean_data.size)

    @property
    def nelem(self) -> int:
        return int(self.quads.shape[0])

    def assemble_data_from_element_coefficients(self, coefficients: FloatArray) -> FloatArray:
        coefficients = np.asarray(coefficients, dtype=np.float64).reshape(-1)
        if coefficients.size != self.nelem:
            raise ValueError(f"Expected {self.nelem} element coefficients")
        weights = self.contribution_base_value * coefficients[self.contribution_element]
        return np.bincount(
            self.assembly_slot,
            weights=weights,
            minlength=self.nnz,
        ).astype(np.float64, copy=False)

    def matrix_from_data(self, data: FloatArray) -> csr_matrix:
        return csr_matrix(
            (np.asarray(data, dtype=np.float64), self.csr_indices, self.csr_indptr),
            shape=(self.n_free, self.n_free),
            copy=False,
        )

    def solve_theta_samples(
        self,
        theta: FloatArray,
        nodal_modal_coefficients: FloatArray,
        *,
        elastic_mean: float,
        elastic_std: float,
        rtol: float = 1e-10,
        atol: float = 1e-14,
        maxiter: int = 80,
        direct_validation_samples: int = 0,
        batch_size: int = 128,
    ) -> tuple[FloatArray, FEMSolveAudit]:
        """Solve all centre deflections for one outer curve.

        The sample stiffness is assembled from the affine KL expansion,
        ``K(theta)=K_mean+sum_j theta_j K_j``.  Each symmetric positive
        definite system is solved by conjugate gradients, preconditioned with
        the exact sparse LU factorisation of the mean stiffness.  Samples that
        do not satisfy the residual check are recomputed by a direct sparse
        solve.  ``batch_size`` is accepted for configuration compatibility;
        the per-sample SciPy kernel is faster for this 1201-DOF model.
        """
        del batch_size
        theta = np.asarray(theta, dtype=np.float64)
        b_nodes = np.asarray(nodal_modal_coefficients, dtype=np.float64)
        if theta.ndim != 2:
            raise ValueError("theta must be n_samples x n_modes")
        if b_nodes.shape != (self.nodes.shape[0], theta.shape[1]):
            raise ValueError(
                "nodal_modal_coefficients must be n_nodes x n_modes and match theta"
            )

        b_elem = b_nodes[self.quads, :].mean(axis=1)
        mode_data = np.empty((theta.shape[1], self.nnz), dtype=np.float64)
        for j in range(theta.shape[1]):
            mode_data[j] = self.assemble_data_from_element_coefficients(
                elastic_std * b_elem[:, j]
            )

        expected_mean = self.assemble_data_from_element_coefficients(
            np.full(self.nelem, elastic_mean, dtype=np.float64)
        )
        if not np.allclose(expected_mean, self.mean_data, rtol=5e-13, atol=1e-8):
            raise RuntimeError("Mean stiffness cache is inconsistent with elastic_mean")

        n_samples = theta.shape[0]
        responses = np.empty(n_samples, dtype=np.float64)
        iteration_counts = np.zeros(n_samples, dtype=np.int32)
        residual_values = np.zeros(n_samples, dtype=np.float64)
        fallback_count = 0
        direct_errors: list[float] = []
        force_norm = max(float(np.linalg.norm(self.force)), np.finfo(np.float64).tiny)
        tiny = np.finfo(np.float64).tiny
        preconditioner = LinearOperator(
            (self.n_free, self.n_free),
            matvec=self._mean_lu.solve,
            dtype=np.float64,
        )

        for sample in range(n_samples):
            data = self.mean_data + theta[sample] @ mode_data
            matrix = self.matrix_from_data(data)
            iterations = 0

            def _count_iteration(_: FloatArray) -> None:
                nonlocal iterations
                iterations += 1

            solution, info = cg(
                matrix,
                self.force,
                x0=self.mean_solution,
                M=preconditioner,
                rtol=rtol,
                atol=atol,
                maxiter=maxiter,
                callback=_count_iteration,
            )
            relative_residual = float(
                np.linalg.norm(matrix @ solution - self.force) / force_norm
            )
            if (
                info != 0
                or not np.all(np.isfinite(solution))
                or not np.isfinite(relative_residual)
                or relative_residual > max(10.0 * rtol, 1e-8)
            ):
                solution = spsolve(matrix.tocsc(), self.force)
                fallback_count += 1
                relative_residual = float(
                    np.linalg.norm(matrix @ solution - self.force) / force_norm
                )

            response = float(self.response_vector @ solution)
            responses[sample] = response
            iteration_counts[sample] = iterations
            residual_values[sample] = relative_residual

            if sample < direct_validation_samples:
                direct = spsolve(matrix.tocsc(), self.force)
                direct_response = float(self.response_vector @ direct)
                direct_errors.append(
                    abs(response - direct_response) / max(abs(direct_response), tiny)
                )

        audit = FEMSolveAudit(
            mean_iterations=float(np.mean(iteration_counts)) if n_samples else 0.0,
            max_iterations=int(np.max(iteration_counts)) if n_samples else 0,
            fallback_count=fallback_count,
            max_relative_residual=float(np.max(residual_values)) if n_samples else 0.0,
            direct_validation_max_relative_error=max(direct_errors, default=0.0),
        )
        return responses, audit


def _generate_mesh(
    a: float, b: float, alpha_deg: float, nx: int, ny: int
) -> tuple[FloatArray, IntArray, IntArray, IntArray, FloatArray, IntArray]:
    alpha = np.deg2rad(alpha_deg)
    nodes = np.empty(((nx + 1) * (ny + 1), 2), dtype=np.float64)
    cursor = 0
    for j in range(ny + 1):
        eta = j / ny
        for i in range(nx + 1):
            xi = i / nx
            nodes[cursor] = [a * xi + b * eta * np.cos(alpha), b * eta * np.sin(alpha)]
            cursor += 1

    quads = np.empty((nx * ny, 4), dtype=np.int64)
    cursor = 0
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            n1 = (j - 1) * (nx + 1) + (i - 1)
            n2 = n1 + 1
            n4 = j * (nx + 1) + (i - 1)
            n3 = n4 + 1
            quads[cursor] = [n1, n2, n3, n4]
            cursor += 1

    all_edges = np.empty((4 * quads.shape[0], 2), dtype=np.int64)
    for e, q in enumerate(quads):
        row = 4 * e
        all_edges[row] = np.sort(q[[0, 1]])
        all_edges[row + 1] = np.sort(q[[1, 2]])
        all_edges[row + 2] = np.sort(q[[2, 3]])
        all_edges[row + 3] = np.sort(q[[3, 0]])
    edges, inverse = np.unique(all_edges, axis=0, return_inverse=True)
    quad_to_edge = inverse.reshape(quads.shape[0], 4)

    edge_normals = np.empty((edges.shape[0], 2), dtype=np.float64)
    for k, (i, j) in enumerate(edges):
        dx = nodes[j, 0] - nodes[i, 0]
        dy = nodes[j, 1] - nodes[i, 1]
        edge_normals[k] = [dy, -dx] / np.hypot(dx, dy)

    bottom = np.arange(nx + 1, dtype=np.int64)
    top = ny * (nx + 1) + np.arange(nx + 1, dtype=np.int64)
    left = np.arange(0, ny * (nx + 1) + 1, nx + 1, dtype=np.int64)
    right = np.arange(nx, (ny + 1) * (nx + 1), nx + 1, dtype=np.int64)
    boundary_nodes = np.unique(np.concatenate((bottom, top, left, right)))
    return nodes, quads, edges, quad_to_edge, edge_normals, boundary_nodes


def _element_data(
    xe: FloatArray,
    ye: FloatArray,
    edge_ids: IntArray,
    edge_normals: FloatArray,
    dbase: FloatArray,
    unit_flex: float,
    q0: float,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    matrix = np.empty((8, 8), dtype=np.float64)
    for i in range(4):
        matrix[i] = poly_row(float(xe[i]), float(ye[i]))
    mids = 0.5 * np.column_stack((xe + np.roll(xe, -1), ye + np.roll(ye, -1)))
    for k in range(4):
        normal = edge_normals[edge_ids[k]]
        matrix[4 + k] = dnormal_row(
            float(mids[k, 0]), float(mids[k, 1]), float(normal[0]), float(normal[1])
        )
    cmat = np.linalg.solve(matrix, np.eye(8, dtype=np.float64))

    ke = np.zeros((8, 8), dtype=np.float64)
    fe = np.zeros(8, dtype=np.float64)
    g = 1.0 / np.sqrt(3.0)
    gauss = ((-g, -g), (g, -g), (g, g), (-g, g))
    for xi, eta in gauss:
        xq, yq, det_j = quad_map(xi, eta, xe, ye)
        bmat = bending_b(xq, yq, cmat)
        nrow = poly_row(xq, yq) @ cmat
        ke += (bmat.T @ dbase @ bmat) * det_j * unit_flex
        fe += q0 * nrow * det_j
    return ke, fe, cmat


def build_plate_model(
    *,
    elastic_mean: float = 2.184e11,
    nu: float = 0.30,
    thickness: float = 0.005,
    a: float = 1.0,
    b: float = 1.0,
    alpha_deg: float = 90.0,
    q0: float = -100.0,
    nx: int = 20,
    ny: int = 20,
) -> PlateModel:
    nodes, quads, edges, quad_to_edge, edge_normals, boundary_nodes = _generate_mesh(
        a, b, alpha_deg, nx, ny
    )
    nnode = nodes.shape[0]
    ndof = nnode + edges.shape[0]
    nelem = quads.shape[0]
    unit_flex = thickness**3 / (12.0 * (1.0 - nu**2))
    dbase = np.array(
        [[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, (1.0 - nu) / 2.0]],
        dtype=np.float64,
    )

    local_ke = np.empty((nelem, 8, 8), dtype=np.float64)
    element_dofs = np.empty((nelem, 8), dtype=np.int64)
    force_full = np.zeros(ndof, dtype=np.float64)
    center_element = int(np.ceil(ny / 2.0) * nx - np.floor(nx / 2.0) - 1)
    center_cmat: FloatArray | None = None

    for e in range(nelem):
        vertices = quads[e]
        edge_ids = quad_to_edge[e]
        xe = nodes[vertices, 0]
        ye = nodes[vertices, 1]
        ke, fe, cmat = _element_data(
            xe, ye, edge_ids, edge_normals, dbase, unit_flex, q0
        )
        dofs = np.concatenate((vertices, nnode + edge_ids))
        local_ke[e] = ke
        element_dofs[e] = dofs
        np.add.at(force_full, dofs, fe)
        if e == center_element:
            center_cmat = cmat

    fixed = boundary_nodes
    free = np.setdiff1d(np.arange(ndof, dtype=np.int64), fixed, assume_unique=True)
    global_to_free = np.full(ndof, -1, dtype=np.int64)
    global_to_free[free] = np.arange(free.size, dtype=np.int64)

    contribution_rows: list[int] = []
    contribution_cols: list[int] = []
    contribution_values: list[float] = []
    contribution_elements: list[int] = []
    for e in range(nelem):
        local_free = global_to_free[element_dofs[e]]
        for a_local in range(8):
            row = int(local_free[a_local])
            if row < 0:
                continue
            for b_local in range(8):
                col = int(local_free[b_local])
                if col < 0:
                    continue
                contribution_rows.append(row)
                contribution_cols.append(col)
                contribution_values.append(float(local_ke[e, a_local, b_local]))
                contribution_elements.append(e)

    rows = np.asarray(contribution_rows, dtype=np.int64)
    cols = np.asarray(contribution_cols, dtype=np.int64)
    base_values = np.asarray(contribution_values, dtype=np.float64)
    element_index = np.asarray(contribution_elements, dtype=np.int64)
    keys = rows * free.size + cols
    unique_keys, slot = np.unique(keys, return_inverse=True)
    csr_rows = unique_keys // free.size
    csr_cols = unique_keys % free.size
    row_counts = np.bincount(csr_rows, minlength=free.size)
    indptr = np.concatenate(([0], np.cumsum(row_counts))).astype(np.int64)
    indices = csr_cols.astype(np.int64, copy=False)

    mean_weights = base_values * elastic_mean
    mean_data = np.bincount(slot, weights=mean_weights, minlength=unique_keys.size)
    mean_matrix = csr_matrix((mean_data, indices, indptr), shape=(free.size, free.size))
    # Enforce the expected symmetry to roundoff without changing the method.
    symmetry_error = np.linalg.norm((mean_matrix - mean_matrix.T).data)
    if symmetry_error > 1e-6 * max(np.linalg.norm(mean_matrix.data), 1.0):
        raise RuntimeError(f"Assembled mean stiffness is unexpectedly nonsymmetric: {symmetry_error}")

    force = force_full[free]
    mean_lu = splu(mean_matrix.tocsc())
    mean_solution = mean_lu.solve(force)

    if center_cmat is None:
        raise RuntimeError("Center element was not initialized")
    alpha = np.deg2rad(alpha_deg)
    xc = 0.5 * a + 0.5 * b * np.cos(alpha)
    yc = 0.5 * b * np.sin(alpha)
    center_vertices = quads[center_element]
    center_edges = quad_to_edge[center_element]
    center_dofs = np.concatenate((center_vertices, nnode + center_edges))
    center_shape = poly_row(float(xc), float(yc)) @ center_cmat
    response_vector = np.zeros(free.size, dtype=np.float64)
    for coefficient, dof in zip(center_shape, center_dofs, strict=True):
        free_id = global_to_free[dof]
        if free_id >= 0:
            response_vector[free_id] += coefficient

    return PlateModel(
        nodes=nodes,
        quads=quads,
        n_free=free.size,
        free_dofs=free,
        force=force,
        response_vector=response_vector,
        csr_indices=indices,
        csr_indptr=indptr,
        assembly_slot=slot.astype(np.int64, copy=False),
        contribution_element=element_index,
        contribution_base_value=base_values,
        mean_data=mean_data,
        mean_matrix=mean_matrix,
        mean_solution=mean_solution,
        _mean_lu=mean_lu,
    )
