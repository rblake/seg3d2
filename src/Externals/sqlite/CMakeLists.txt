cmake_minimum_required(VERSION 2.6)
PROJECT(SQLite)


###########################################
# Source files for SQLite
SET( SQLITE_SRCS
	sqlite3.h
	sqlite3ext.h
	shell.c
	sqlite3.c
)

###########################################
# Set the name of the library

SET(SCI_SQLITE_LIBRARY sqlite CACHE INTERNAL "SQLite Library name." FORCE)

###########################################
# Include directory of SQLite
SET(SCI_SQLITE_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "SQLite include directories" FORCE)

ADD_LIBRARY(${SCI_SQLITE_LIBRARY} STATIC ${SQLITE_SRCS})