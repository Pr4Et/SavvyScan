SavvyscanPlugin.dll should be added to the SerialEM installation folder to support connection with the server.
The new plugin is decoupled from SerialEM resources and is thus suitbale for all SerialEM versions.
The installtion requires definition of the server ip address outside the SeriaEMproperties.txt file as follows:
An environment variable SAVVY_SERVER_IP should be set with the server IP address (In Windows: go to Advanced System Settings, Environment Varaibles). The default address is 192.168.100.90.
 