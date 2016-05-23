#!/bin/sh
case `uname` in Darwin*) glibtoolize --copy ;;
  *) libtoolize --copy ;; esac
autoreconf --force --install -v
