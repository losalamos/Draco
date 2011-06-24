macro( draco_flat_function_names )

   if( ${SYSTEM} MATCHES "ibm-aix" )
       set(DACS_DEVICE_INIT dacs_device_init)
       set(DACS_DEVICE_GET_DE_ID dacs_device_get_de_id)
       set(DACS_DEVICE_GET_PID dacs_device_get_pid)

   else()
       set(DACS_DEVICE_INIT dacs_device_init_)
       set(DACS_DEVICE_GET_DE_ID dacs_device_get_de_id_)
       set(DACS_DEVICE_GET_PID dacs_device_get_pid_)

   endif()

endmacro()
