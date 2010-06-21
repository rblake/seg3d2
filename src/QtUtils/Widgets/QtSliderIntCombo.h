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

#ifndef QTUTILS_WIDGETS_QTSLIDERINTCOMBO_H
#define QTUTILS_WIDGETS_QTSLIDERINTCOMBO_H

// QT Includes
#include <QWidget>

// Boost includes
#include <boost/shared_ptr.hpp>

namespace QtUtils 
{

// Forward declaration
class QtSliderIntCombo;
class QtSliderIntComboPrivate;
typedef boost::shared_ptr< QtSliderIntComboPrivate > QtSliderIntComboPrivateHandle;

// Class definition
class QtSliderIntCombo : public QWidget
{
Q_OBJECT

// -- constructor/destructor --
public:
    QtSliderIntCombo( QWidget* parent = 0, bool edit_range = false );
    virtual ~QtSliderIntCombo();

Q_SIGNALS:
	void valueAdjusted( int );
	void rangeChanged( int, int );
    
public Q_SLOTS:
    void setStep(int);
	void setRange( int, int );
	void setCurrentValue( int );
	
public:
	int get_value() { return value_; }
//	void block() { blockSignals( true ); }
//	void unblock() { blockSignals( false ); }
//	void parent(){ return parentWidget(); }
	

// -- widget internals -- 
private:
    QtSliderIntComboPrivateHandle private_;
    
private Q_SLOTS:
    void edit_ranges( bool edit );
    void change_min( int new_min );
    void change_max( int new_max );
    void double_range();
    void half_range();
    void slider_signal( int value );
    void spinner_signal( int value );
    
private:
	void block_signals( bool block );
	int value_;
  
};

}  // end namespace QtUtils

#endif
