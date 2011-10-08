#!/bin/bash

all: rmap cgpost

rmap:
	cd reverse-map/src && $(MAKE)

cgpost:
	cd cg-post/src && $(MAKE)

clean:
	cd cg-post/src && $(MAKE) clean
	cd reverse-map/src && $(MAKE) clean
