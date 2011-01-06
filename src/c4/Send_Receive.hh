//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Send_Receive.hh
 * \author Mike Buksas
 * \brief  Define class Send_Receive
 * \note   Copyright (C) 2007-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_Send_Receive_hh
#define c4_Send_Receive_hh

#include "C4_Functions.hh"
#include "C4_Req.hh"
#include "ds++/Assert.hh"

namespace rtt_c4
{

//===========================================================================//
/*!
 * \class Utitlity classes for variable-size communication of data
 *
 * \brief These classes facilitate the construction of special-purpose data
 * communicators between processors by encapsulating an asynchnrous send and
 * synchronous receive of the size and data. The operations are templated on
 * the data type.
 *
 * To use these classes, derive from them and use the protected member
 * functions. The derived Sender class can generate or accquire the
 * information to send in any desired way (e.g. extracting it out of a larger
 * set of data) and use the Sender::send function to send the size and then
 * the actual data.
 *
 * A derived class of Receiver has two methods for receiving data: It can get
 * the size (Receiver::receive_size) and allocate it's own storage for the
 * data (Receiver::receive_data), or use Receiver::recieve which allocates
 * array storage for the data and returns the size. The allocated data must be
 * managed by the calling code. (Generally, the derived class)
 *
 * The base class send/receive functions obey the convention that data of size
 * zero is not actually sent.
 *
 */
/*! 
 * \example c4/test/tstSend_Receive.cc
 *
 * Test of Send_Receive.
 */
//===========================================================================//

const int SIZE_CHANNEL = 325;
const int DATA_CHANNEL = 326;

class Sender
{
  private:

    int to_node;
    C4_Req size_handle, data_handle;
    int sz;

    
  public:
  
    explicit Sender(int node) 
        : to_node(node),
          size_handle( C4_Req() ),
          data_handle( C4_Req() )
    {
        Check (to_node >= 0);
        Check (to_node < rtt_c4::nodes());
    }
    virtual ~Sender() {/* empty */}
    
  protected:

    template <typename T>
    void send(int size, T* data)
    {
        
        // Hold onto size in a private member variable, to avoid
        // having size fall out of scope before the send actually
        // completes.
        sz = size;
        size_handle = send_async(&sz, 1, to_node, SIZE_CHANNEL);
        if (size > 0)
        {
            data_handle = rtt_c4::send_async(data, size, to_node, DATA_CHANNEL);
        }

    }

    void wait()
    {
        size_handle.wait();
        data_handle.wait();
    }

};



class Receiver
{

  private:

    int from_node;

  public:

    explicit Receiver(int node) 
      : from_node(node)
    {
        Check (from_node >= 0);
        Check (from_node < rtt_c4::nodes());
    }
    virtual ~Receiver() {/*empty*/}

  protected:

    int receive_size() const
    {        
        int size_rec(0);
        /*int count = */ rtt_c4::receive<int>(&size_rec, 1, from_node, SIZE_CHANNEL);

        return size_rec;
    }


    template <typename T>
    void receive_data(int size, T* data) const
    {
        if (size > 0)
        {
            Check(data);
            rtt_c4::receive<T>(data, size, from_node, DATA_CHANNEL);

        }
    }

    
    template <typename T>
    int receive(T* &data) const
    {
        int size = receive_size();
        if (size > 0)
        {
            data = new T[size];
            receive_data(size, data);
        }

        return size;
    }

};

} // end namespace rtt_c4

#endif // c4_Send_Receive_hh

//---------------------------------------------------------------------------//
//              end of c4/Send_Receive.hh
//---------------------------------------------------------------------------//
