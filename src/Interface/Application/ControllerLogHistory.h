/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2009 Scientific Computing and Imaging Institute,
 University of Utah.
 
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
 */

#ifndef INTERFACE_APPLICATION_CONTROLLERLOGHISTORY_H
#define INTERFACE_APPLICATION_CONTROLLERLOGHISTORY_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

#ifndef Q_MOC_RUN

// STL includes
#include <string>
#include <deque>

// QT includes
#include <QtCore/QObject>
#include <QtCore/QVariant>
#include <QtCore/QModelIndex>

// Core includes
#include <Core/Utils/Log.h>

#endif

namespace Seg3D
{

class ControllerLogHistory : public QAbstractTableModel
{

Q_OBJECT

public:
	ControllerLogHistory( size_t log_size, QObject* parent = 0 );

	virtual ~ControllerLogHistory();

	int rowCount( const QModelIndex &index ) const;
	int columnCount( const QModelIndex &index ) const;

	QVariant data( const QModelIndex& index, int role ) const;
	QVariant headerData( int section, Qt::Orientation orientation, int role ) const;

	void add_log_entry( int message_type, std::string& message );

private:
	// Classes needed for storing the recent log history
	typedef std::pair< int, std::string > log_entry_type;
	typedef std::deque< log_entry_type > log_history_type;

	log_history_type log_history_;
	size_t log_history_size_;
};

} // end namespace Seg3D

#endif
