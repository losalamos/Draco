//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde_pt.cc
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  Explicit template instatiations for class kde.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved.
 *
 */
//------------------------------------------------------------------------------------------------//

#include "kde.t.hh"

namespace kde {

// Explicit template instantiations go here.
template class kde<kde_coordinates::CART>;

} // namespace  kde

//------------------------------------------------------------------------------------------------//
// end of kde_pt.cc
//------------------------------------------------------------------------------------------------//
