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

#ifndef QTUTILS_BRIDGE_DETAIL_QTTRANSFERFUNCTIONCURVECONNECTOR_H
#define QTUTILS_BRIDGE_DETAIL_QTTRANSFERFUNCTIONCURVECONNECTOR_H

#include <QPointer>

#include <Core/VolumeRenderer/TransferFunctionFeature.h>

#include <QtUtils/Widgets/QtTransferFunctionCurve.h>
#include <QtUtils/Bridge/detail/QtConnectorBase.h>

namespace QtUtils
{

class QtTransferFunctionCurveConnector : public QtConnectorBase
{
	Q_OBJECT

public:
	QtTransferFunctionCurveConnector( QtTransferFunctionCurve* parent, 
		Core::TransferFunctionFeatureHandle& tf_feature, bool blocking = true );
	virtual ~QtTransferFunctionCurveConnector();

	// -- slot functions for boost signals --
private:
	static void SetCurveControlPoints( QPointer< QtTransferFunctionCurveConnector > qpointer,
		Core::TransferFunctionControlPointVector control_points, Core::ActionSource source );

	static void UpdateCurveColor( QPointer< QtTransferFunctionCurveConnector > qpointer );

	// -- slot functions for Qt signals --
private Q_SLOTS:
	void set_control_points_state( const Core::TransferFunctionControlPointVector& control_points );

private:
	QtTransferFunctionCurve* parent_;
	Core::TransferFunctionFeatureHandle tf_feature_;
};

}

#endif