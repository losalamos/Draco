//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Surface_Descriptor.cc
 * \author Mike Buksas
 * \date   Tue Aug 12 15:17:31 2003
 * \brief  Implementation file for Surface_Descriptor
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Surface_Descriptor.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"

using namespace std;

namespace rtt_mc
{

//---------------------------------------------------------------------------//
const int Surface_Descriptor::sizes[1] = {2};

//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor
 * 
 * \param int type of the surface
 * \param data for the surface
 */
Surface_Descriptor::Surface_Descriptor(int type_, const std::vector<double>& data_)
    : type(static_cast<Surface_Descriptor::Surface_Type>(type_)), data(data_)
{

    Check(type >= 0); Check(type < kinds);

    Check(data.size() == sizes[type]);

}


//---------------------------------------------------------------------------//
/*! 
 * \brief Unpcaking consutrctor
 * 
 * \param data vector of char data
 */
Surface_Descriptor::Surface_Descriptor(const std::vector<char>& packed)
{

    rtt_dsxx::Unpacker u;

    u.set_buffer(packed.size(), &packed[0]);

    u >> type;  Require(type >= 0); Require(type < kinds);

    data.resize( sizes[type] );

    for (int i = 0; i < sizes[type]; ++i) { u >> data[i]; }

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Packing Operator
 * 
 * \return a vector<char> containing the data of the object
 */
vector<char> Surface_Descriptor::pack() const
{

    rtt_dsxx::Packer p;

    int size = sizeof(int) + sizeof(double) * data.size();

    vector<char> packed(size);

    p.set_buffer(size, &packed[0]);

    p << type;

    for (int i=0; i<data.size(); ++i) p << data[i];

    return packed;

}

//---------------------------------------------------------------------------//
bool Surface_Descriptor::operator==(const Surface_Descriptor& rhs) const
{
    return (data == rhs.data && type == rhs.type);
}

bool Surface_Descriptor::operator!=(const Surface_Descriptor& rhs) const
{
    return !(*this == rhs);
}


} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                 end of Surface_Descriptor.cc
//---------------------------------------------------------------------------//
