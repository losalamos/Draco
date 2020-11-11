//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.t.hh
 * \author Mathew Cleveland
 * \date   November 10th 2020
 * \brief  Templated KDE functions for various dimensions and coordinate
 * systems
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef kde_kde_t_hh
#define kde_kde_t_hh

#include "kde.hh"
#include <iostream>

namespace kde {

/*!
 * reconstruction 
 * \brief
 *
 * Cartesian geometry reconstruction of a 1D distribution
 *
 * \param[in] distribution 
 * \param[in] position
 * \param[in] band_width
 * \param[inout] final_distribution returned final function distribution
 *
 * Test of kde.
 */ 
template<>
template<>
void kde<kde_coordinates::CART>::reconstruction<1>(const std::vector<double> &distribution, const std::vector<double> position, const std::vector<double> band_width, std::vector<double> &final_distribution) const{
    Check(position.size()==distribution.size());
    Check(band_width.size()==distribution.size());
    Check(final_distribution.size()==distribution.size());
    std::cout<<"HELLO WORLD"<<std::endl;
    std::cout<<"Coordinate = CART"<<std::endl;
    std::cout<<"DIM = 1"<<std::endl;
}

/*!
 * reconstruction 
 * \brief
 *
 * Cartesian geometry reconstruction of a 2D distribution
 *
 * \param[in] distribution 
 * \param[in] position
 * \param[in] band_width
 * \param[inout] final_distribution returned final function distribution
 *
 * Test of kde.
 */ 
template<>
template<>
void kde<kde_coordinates::CART>::reconstruction<2>(const std::vector<double> &distribution, const std::vector<double> position, const std::vector<double> band_width, std::vector<double> &final_distribution) const{
    Check(position.size()==distribution.size());
    Check(band_width.size()==distribution.size());
    Check(final_distribution.size()==distribution.size());
    std::cout<<"HELLO WORLD"<<std::endl;
    std::cout<<"Coordinate = CART"<<std::endl;
    std::cout<<"DIM = 2"<<std::endl;
}


} // end namespace  kde

#endif // kde_kde_t_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.t.hh
//------------------------------------------------------------------------------------------------//

