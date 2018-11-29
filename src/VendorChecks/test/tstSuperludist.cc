//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   VendorChecks/test/tstSuperludist.cc
 * \date   Monday, May 16, 2016, 16:30 pm
 * \brief  Attempt to link to libsuperludist and run a simple problem.
 * \note   Copyright (C) 2016-2018, Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This code is a modified version of \c pddrive.c provided in the EXAMPLES
 * directory distributed with the SuperLU_DIST-4.3 source. It is the driver
 * program for the PDGSSVX example. The original code is dated 11/1/2007 and is
 * attributed to Lawrence Berkeley National Lab, Univ. of California Berkeley.
 */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/path.hh"
#include <sstream>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundef"
#include <superlu_ddefs.h>
#pragma GCC diagnostic pop
#include <vector>

// forward declarations
void test_superludist(rtt_c4::ParallelUnitTest &ut);
int dcreate_matrix(SuperMatrix *A, int nrhs, double **rhs, int *ldb, double **x,
                   int *ldx, FILE *fp, gridinfo_t *grid);

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_superludist(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
/* Purpose
 * =======
 *
 * The driver program PDDRIVE.
 *
 * This example illustrates how to use PDGSSVX with the full (default) options
 * to solve a linear system.
 *
 * Five basic steps are required:
 *   1. Initialize the MPI environment and the SuperLU process grid
 *   2. Set up the input matrix and the right-hand side
 *   3. Set the options argument
 *   4. Call pdgssvx
 *   5. Release the process grid and terminate the MPI environment
 */
void test_superludist(rtt_c4::ParallelUnitTest &ut) {
// for superlu-dist > version 4, SUPERLU_DIST_MAJOR_VERSION is defined in
// superlu_defs.h
#ifdef SUPERLU_DIST_MAJOR_VERSION
  superlu_dist_options_t options;
#else
  superlu_options_t options;
#endif
  SuperLUStat_t stat;
  SuperMatrix A;
  ScalePermstruct_t ScalePermstruct;
  LUstruct_t LUstruct;
  SOLVEstruct_t SOLVEstruct;
  gridinfo_t grid;
  double *berr;
  double *b;
  double *xtrue;
  int m(-42);
  int n(-42);
  int const nprow = 2; // Must run mpi with nprow*npcol ranks.
  int const npcol = 2;
  int iam(-42);
  int info(-42);
  int ldb(-42);
  int ldx(-42);
  int nrhs(1); // Number of right-hand side.
  std::string const inpPath = ut.getTestSourcePath();
  std::string const filename(inpPath + "big.rua");
  FILE *fp;

  //---------------------------------------------------------------------------//
  // Initialize the SuperLU process grid
  //---------------------------------------------------------------------------//

  // This file contains the matrix (g20.rua)
  fp = fopen(filename.c_str(), "r");
  Insist(fp != NULL, "File g20.rua does not exist or cannot be loaded.");

  superlu_gridinit(rtt_c4::communicator, nprow, npcol, &grid);

  // Bail out if I do not belong in the grid.
  iam = grid.iam;
  if (iam >= nprow * npcol) {
    superlu_gridexit(&grid);
    return;
  }

  if (!iam) {
    std::cout << "\nInput matrix file: " << filename;
    std::cout << "\nProcess grid: " << static_cast<int>(grid.nprow) << " X "
              << static_cast<int>(grid.npcol) << "\n"
              << std::endl;
  }

  //---------------------------------------------------------------------------//
  // Load the matirx from file and setup the RHS.
  //---------------------------------------------------------------------------//

  dcreate_matrix(&A, nrhs, &b, &ldb, &xtrue, &ldx, fp, &grid);

  berr = doubleMalloc_dist(nrhs);
  Insist(berr != NULL, "Malloc fails for berr[].");

  if (A.Stype != SLU_NR_loc)
    ITFAILS;
  if (A.Dtype != SLU_D)
    ITFAILS;
  if (A.Mtype != SLU_GE)
    ITFAILS;
  if (A.nrow != 4960)
    ITFAILS;
  if (A.ncol != 4960)
    ITFAILS;
  PASSMSG("Successfully loaded data file.");

  //---------------------------------------------------------------------------//
  // Solve the linear system
  //---------------------------------------------------------------------------//

  /* Set the default input options:
       options.Fact              = DOFACT;
       options.Equil             = YES;
       options.ParSymbFact       = NO;
       options.ColPerm           = METIS_AT_PLUS_A;
       options.RowPerm           = LargeDiag;      // version < 5.4
       options.RowPerm           = LargeDiag_MC64; // version > 5.4
       options.ReplaceTinyPivot  = YES;
       options.IterRefine        = DOUBLE;
       options.Trans             = NOTRANS;
       options.SolveInitialized  = NO;
       options.RefineInitialized = NO;
       options.PrintStat         = YES;
    */
  set_default_options_dist(&options);

  if (options.Fact != DOFACT)
    ITFAILS;
  if (options.Equil != YES)
    ITFAILS;
  if (options.ColPerm != METIS_AT_PLUS_A)
    ITFAILS;
  if (options.IterRefine != SLU_DOUBLE)
    ITFAILS;
#if SUPERLU_DIST_MAJOR_VERSION == 5 && SUPERLU_DIST_MINOR_VERSION >= 4
  if (options.RowPerm != LargeDiag_MC64)
    ITFAILS;
#else
  if (options.RowPerm != LargeDiag)
    ITFAILS;
#endif

  if (!iam) {
    print_sp_ienv_dist(&options);
    print_options_dist(&options);
  }

  m = A.nrow;
  n = A.ncol;

  // Initialize ScalePermstruct and LUstruct.
  ScalePermstructInit(m, n, &ScalePermstruct);
  LUstructInit(n, &LUstruct);

  if (ScalePermstruct.DiagScale != NOEQUIL)
    ITFAILS;

  // Initialize the statistics variables.
  PStatInit(&stat);

  if (!rtt_dsxx::soft_equiv(static_cast<double>(*stat.ops), 0.0))
    ITFAILS;
  if (stat.TinyPivots != 0)
    ITFAILS;

  // Call the linear equation solver.
  pdgssvx(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid, &LUstruct,
          &SOLVEstruct, berr, &stat, &info);

  if (!rtt_dsxx::soft_equiv(*berr, 0.0))
    ITFAILS;
  if (info != 0)
    ITFAILS;
  if (options.SolveInitialized != YES)
    ITFAILS;
  if (options.RefineInitialized != YES)
    ITFAILS;
  if (ScalePermstruct.DiagScale != BOTH)
    ITFAILS;

  // Check the accuracy of the solution.
  pdinf_norm_error(iam, ((NRformat_loc *)A.Store)->m_loc, nrhs, b, ldb, xtrue,
                   ldx, &grid);

  /* Print the statistics. */
  PStatPrint(&options, &stat, &grid);

  //---------------------------------------------------------------------------//
  // Deallocate storage
  //---------------------------------------------------------------------------//

  PStatFree(&stat);
  Destroy_CompRowLoc_Matrix_dist(&A);
  ScalePermstructFree(&ScalePermstruct);
  Destroy_LU(n, &grid, &LUstruct);
  LUstructFree(&LUstruct);
  if (options.SolveInitialized)
    dSolveFinalize(&options, &SOLVEstruct);

  SUPERLU_FREE(b);
  SUPERLU_FREE(xtrue);
  SUPERLU_FREE(berr);

  //---------------------------------------------------------------------------//
  // Release the SuperLU process grid
  //---------------------------------------------------------------------------//

  superlu_gridexit(&grid);

  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read the matrix from data file in Harwell-Boeing format, and
 *        distribute it to processors in a distributed compressed row format. It
 *        also generate the distributed true solution X and the right-hand side
 *        RHS.
 *
 * \param[out] A    local matrix A in NR_loc format.
 * \param[in]  nrhs number of right-hand sides.
 * \param[out] rhs  the right-hand side matrix.
 * \param[out] ldb  leading dimension of the right-hand side matrix.
 * \param[out] x    the true solution matrix.
 * \param[out] ldx  the leading dimension of the true solution matrix.
 * \param[in]  fp   the matrix file pointer.
 * \param[in]  grid the 2D process mesh
 * \return 0
 *
 * This function is a modified version of \c dcreate_matrix.c found in the
 * EXAMPLES directory of the SuperLU_DIST source distribution.
 */
int dcreate_matrix(SuperMatrix *A, int nrhs, double **rhs, int *ldb, double **x,
                   int *ldx, FILE *fp, gridinfo_t *grid) {
  SuperMatrix GA;                  /* global A */
  double *b_global, *xtrue_global; /* replicated on all processes */
  int_t *rowind, *colptr;          /* global */
  double *nzval;                   /* global */
  double *nzval_loc;               /* local */
  int_t *colind, *rowptr;          /* local */
  int_t m, n, nnz;
  int_t m_loc, fst_row, nnz_loc;
  int_t m_loc_fst; /* Record m_loc of the first p-1 processors,
                           when mod(m, p) is not zero. */
  int_t row, i, j, relpos;
  int iam;
  char trans[1];
  int_t *marker;

  iam = grid->iam;

  if (!iam) {
    /* Read the matrix stored on disk in Harwell-Boeing format. */
    dreadhb_dist(iam, fp, &m, &n, &nnz, &nzval, &rowind, &colptr);

    /* Broadcast matrix A to the other PEs. */
    MPI_Bcast(&m, 1, mpi_int_t, 0, grid->comm);
    MPI_Bcast(&n, 1, mpi_int_t, 0, grid->comm);
    MPI_Bcast(&nnz, 1, mpi_int_t, 0, grid->comm);
    MPI_Bcast(nzval, nnz, MPI_DOUBLE, 0, grid->comm);
    MPI_Bcast(rowind, nnz, mpi_int_t, 0, grid->comm);
    MPI_Bcast(colptr, n + 1, mpi_int_t, 0, grid->comm);
  } else {
    /* Receive matrix A from PE 0. */
    MPI_Bcast(&m, 1, mpi_int_t, 0, grid->comm);
    MPI_Bcast(&n, 1, mpi_int_t, 0, grid->comm);
    MPI_Bcast(&nnz, 1, mpi_int_t, 0, grid->comm);

    /* Allocate storage for compressed column representation. */
    dallocateA_dist(n, nnz, &nzval, &rowind, &colptr);

    MPI_Bcast(nzval, nnz, MPI_DOUBLE, 0, grid->comm);
    MPI_Bcast(rowind, nnz, mpi_int_t, 0, grid->comm);
    MPI_Bcast(colptr, n + 1, mpi_int_t, 0, grid->comm);
  }

  /* Compute the number of rows to be distributed to local process */
  m_loc = m / (grid->nprow * grid->npcol);
  m_loc_fst = m_loc;
  /* When m / procs is not an integer */
  if ((m_loc * grid->nprow * grid->npcol) != m) {
    if (iam == (grid->nprow * grid->npcol - 1)) /* last proc. gets all*/
      m_loc = m - m_loc * (grid->nprow * grid->npcol - 1);
  }

  /* Create compressed column matrix for GA. */
  dCreate_CompCol_Matrix_dist(&GA, m, n, nnz, nzval, rowind, colptr, SLU_NC,
                              SLU_D, SLU_GE);

  /* Generate the exact solution and compute the right-hand side. */
  if (!(b_global = doubleMalloc_dist(m * nrhs)))
    ABORT("Malloc fails for b[]");
  if (!(xtrue_global = doubleMalloc_dist(n * nrhs)))
    ABORT("Malloc fails for xtrue[]");
  *trans = 'N';

  dGenXtrue_dist(n, nrhs, xtrue_global, n);
  dFillRHS_dist(trans, nrhs, xtrue_global, n, &GA, b_global, m);

  /*************************************************
     * Change GA to a local A with NR_loc format     *
     *************************************************/

  rowptr = (int_t *)intMalloc_dist(m_loc + 1);
  marker = (int_t *)intCalloc_dist(n);

  /* Get counts of each row of GA */
  for (i = 0; i < n; ++i)
    for (j = colptr[i]; j < colptr[i + 1]; ++j)
      ++marker[rowind[j]];
  /* Set up row pointers */
  rowptr[0] = 0;
  fst_row = iam * m_loc_fst;
  nnz_loc = 0;
  for (j = 0; j < m_loc; ++j) {

    row = fst_row + j;
    rowptr[j + 1] = rowptr[j] + marker[row];
    marker[j] = rowptr[j];
  }
  nnz_loc = rowptr[m_loc];

  nzval_loc = (double *)doubleMalloc_dist(nnz_loc);
  colind = (int_t *)intMalloc_dist(nnz_loc);

  /* Transfer the matrix into the compressed row storage */
  for (i = 0; i < n; ++i) {
    for (j = colptr[i]; j < colptr[i + 1]; ++j) {
      row = rowind[j];
      if ((row >= fst_row) && (row < fst_row + m_loc)) {
        row = row - fst_row;
        relpos = marker[row];
        colind[relpos] = i;
        nzval_loc[relpos] = nzval[j];
        ++marker[row];
      }
    }
  }

  /* Destroy GA */
  Destroy_CompCol_Matrix_dist(&GA);

  /******************************************************/
  /* Change GA to a local A with NR_loc format */
  /******************************************************/

  /* Set up the local A in NR_loc format */
  dCreate_CompRowLoc_Matrix_dist(A, m, n, nnz_loc, m_loc, fst_row, nzval_loc,
                                 colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE);

  /* Get the local B */
  if (!((*rhs) = doubleMalloc_dist(m_loc * nrhs)))
    ABORT("Malloc fails for rhs[]");
  for (j = 0; j < nrhs; ++j) {
    for (i = 0; i < m_loc; ++i) {
      row = fst_row + i;
      (*rhs)[j * m_loc + i] = b_global[j * n + row];
    }
  }
  *ldb = m_loc;

  /* Set the true X */
  *ldx = m_loc;
  if (!((*x) = doubleMalloc_dist(*ldx * nrhs)))
    ABORT("Malloc fails for x_loc[]");

  /* Get the local part of xtrue_global */
  for (j = 0; j < nrhs; ++j) {
    for (i = 0; i < m_loc; ++i)
      (*x)[i + j * (*ldx)] = xtrue_global[i + fst_row + j * n];
  }

  SUPERLU_FREE(b_global);
  SUPERLU_FREE(xtrue_global);
  SUPERLU_FREE(marker);

  return 0;
}

//---------------------------------------------------------------------------//
// end of tstSuperludist.cc
//---------------------------------------------------------------------------//
