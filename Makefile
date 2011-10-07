#!/bin/bash

all:
	cd cg-post/src && $(MAKE)

clean:
	cd cg-post/src && $(MAKE) clean
