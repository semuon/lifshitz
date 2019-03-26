#ifndef OUR_BOOST_SMART_PTR_MAKE_SHARED_HPP_INCLUDED
#define OUR_BOOST_SMART_PTR_MAKE_SHARED_HPP_INCLUDED


#include <boost/config.hpp>
#include <boost/move/core.hpp>
#include <boost/move/utility_core.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/detail/sp_forward.hpp>
#include <boost/type_traits/type_with_alignment.hpp>
#include <boost/type_traits/alignment_of.hpp>
#include <cstddef>
#include <new>


#include <boost/smart_ptr/make_shared_object.hpp>



#if !defined( BOOST_NO_FUNCTION_TEMPLATE_ORDERING )
# define BOOST_SP_MSD( T ) boost::detail::sp_inplace_tag< boost::detail::sp_ms_deleter< T > >()
#else
# define BOOST_SP_MSD( T ) boost::detail::sp_ms_deleter< T >()
#endif

template< class T > typename boost::detail::sp_if_not_array< T >::type our_make_shared()
{
    return boost::make_shared<T>();    
}


// C++03 version

template< class T, class A1 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1 )
{
    return boost::make_shared<T>(a1); 
}


template< class T, class A1, class A2 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2 )
{
    return boost::make_shared<T>(a1. a2);
}



template< class T, class A1, class A2, class A3 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3 )
{
   return boost::make_shared<T>(a1, a2, a3);
}



template< class T, class A1, class A2, class A3, class A4 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4 )
{
    return boost::make_shared<T>(a1, a2, a3, a4);
}



template< class T, class A1, class A2, class A3, class A4, class A5 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5 )
{
  return boost::make_shared<T>(a1, a2, a3 ,a4, a5);
}


template< class T, class A1, class A2, class A3, class A4, class A5, class A6 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5, BOOST_FWD_REF(A6) a6 )
{
   return boost::make_shared<T>(a1, a2, a3 ,a4, a5, a6);
}



template< class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5, BOOST_FWD_REF(A6) a6, BOOST_FWD_REF(A7) a7 )
{
    return boost::make_shared<T>(a1, a2, a3 ,a4, a5, a6, a7);
}



template< class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5, BOOST_FWD_REF(A6) a6, BOOST_FWD_REF(A7) a7, BOOST_FWD_REF(A8) a8 )
{
   return boost::make_shared<T>(a1, a2, a3 ,a4, a5, a6, a7, a8);
}



template< class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9 >
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5, BOOST_FWD_REF(A6) a6, BOOST_FWD_REF(A7) a7, BOOST_FWD_REF(A8) a8, BOOST_FWD_REF(A9) a9 )
{
    return boost::make_shared<T>(a1, a2, a3 ,a4, a5, a6, a7, a8 , a9);
}



template< class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12 , class A13, class A14 , class A15>
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5, BOOST_FWD_REF(A6) a6, BOOST_FWD_REF(A7) a7, BOOST_FWD_REF(A8) a8, BOOST_FWD_REF(A9) a9,  BOOST_FWD_REF(A10) a10, BOOST_FWD_REF(A11) a11, BOOST_FWD_REF(A12) a12,  BOOST_FWD_REF(A13) a13, BOOST_FWD_REF(A14) a14, BOOST_FWD_REF(A15) a15  )
{
    boost::shared_ptr< T > pt( static_cast< T* >( 0 ), BOOST_SP_MSD( T ) );

    boost::detail::sp_ms_deleter< T > * pd = static_cast<boost::detail::sp_ms_deleter< T > *>( pt._internal_get_untyped_deleter() );

    void * pv = pd->address();

    ::new( pv ) T(
        boost::forward<A1>( a1 ),
        boost::forward<A2>( a2 ),
        boost::forward<A3>( a3 ),
        boost::forward<A4>( a4 ),
        boost::forward<A5>( a5 ),
        boost::forward<A6>( a6 ),
        boost::forward<A7>( a7 ),
        boost::forward<A8>( a8 ),
        boost::forward<A9>( a9 ),
        boost::forward<A10>( a10 ),
	boost::forward<A11>( a11 ),
	boost::forward<A12>( a12 ),
	boost::forward<A13>( a13 ),
	boost::forward<A14>( a14 ),
        boost::forward<A15>( a15 )
		  
        );

    pd->set_initialized();

    T * pt2 = static_cast< T* >( pv );

    boost::detail::sp_enable_shared_from_this( &pt, pt2, pt2 );
    return boost::shared_ptr< T >( pt, pt2 );
}



template< class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11, class A12 , class A13, class A14, class A15, class A16>
typename boost::detail::sp_if_not_array< T >::type our_make_shared( BOOST_FWD_REF(A1) a1, BOOST_FWD_REF(A2) a2, BOOST_FWD_REF(A3) a3, BOOST_FWD_REF(A4) a4, BOOST_FWD_REF(A5) a5, BOOST_FWD_REF(A6) a6, BOOST_FWD_REF(A7) a7, BOOST_FWD_REF(A8) a8, BOOST_FWD_REF(A9) a9,  BOOST_FWD_REF(A10) a10, BOOST_FWD_REF(A11) a11, BOOST_FWD_REF(A12) a12,  BOOST_FWD_REF(A13) a13, BOOST_FWD_REF(A14) a14, BOOST_FWD_REF(A15) a15, BOOST_FWD_REF(A16) a16 )
{
    boost::shared_ptr< T > pt( static_cast< T* >( 0 ), BOOST_SP_MSD( T ) );

    boost::detail::sp_ms_deleter< T > * pd = static_cast<boost::detail::sp_ms_deleter< T > *>( pt._internal_get_untyped_deleter() );

    void * pv = pd->address();

    ::new( pv ) T(
        boost::forward<A1>( a1 ),
        boost::forward<A2>( a2 ),
        boost::forward<A3>( a3 ),
        boost::forward<A4>( a4 ),
        boost::forward<A5>( a5 ),
        boost::forward<A6>( a6 ),
        boost::forward<A7>( a7 ),
        boost::forward<A8>( a8 ),
        boost::forward<A9>( a9 ),
        boost::forward<A10>( a10 ),
	boost::forward<A11>( a11 ),
	boost::forward<A12>( a12 ),
	boost::forward<A13>( a13 ),
	boost::forward<A14>( a14 ), 
	boost::forward<A15>( a15 ),
	boost::forward<A16>( a16 )
		  
        );

    pd->set_initialized();

    T * pt2 = static_cast< T* >( pv );

    boost::detail::sp_enable_shared_from_this( &pt, pt2, pt2 );
    return boost::shared_ptr< T >( pt, pt2 );
}



#endif 