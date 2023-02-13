This directory contains examples for the Embedded Server Option.

Client
Server
These two programs implement a very simple client-server communication. The server is intended to run on the Netbox while the client runs on the user's PC.

dbus
This is an example on how to connect to the Netbox's internal signals (currently only LAN state).

simple_rec_fifo_mail
This example will run a FIFO multi acquisition and send a mail for each acquired segment as a SBench6-compatible binary file. Please keep in mind that a high trigger frequency will flood your mailserver with emails which might trigger some spam detection mechanisms.

