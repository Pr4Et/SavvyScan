#
# Define a list of pkg-config packages we want to use
pkg_packages := dbus-glib-1
PKG_CFLAGS  := $(shell pkg-config --cflags $(pkg_packages))
PKG_LDFLAGS := $(shell pkg-config --libs $(pkg_packages))
# Additional flags for the compiler:
#    -g : Add debugging symbols
# -Wall : Enable most gcc warnings
ADD_CFLAGS := -g -Wall
# Combine user supplied, additional, and pkg-config flags
CFLAGS  := $(PKG_CFLAGS) $(ADD_CFLAGS) $(CFLAGS)
LDFLAGS := $(PKG_LDFLAGS) $(LDFLAGS)
targets = dbus-example
.PHONY: all clean
all: $(targets)
dbus-example: dbus-example.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
clean:
	$(RM) $(targets)
