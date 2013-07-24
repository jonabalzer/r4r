/*////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the R4R library.
//
// The R4R library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The R4R library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////*/

#ifndef PREFERENCES_H
#define PREFERENCES_H

#include <QWidget>

#include "params.h"

using namespace R4R;
using namespace std;

namespace Ui {
class Preferences;
}

class Preferences : public QWidget
{
    Q_OBJECT
    
public:
    explicit Preferences(QWidget *parent = 0);
    ~Preferences();
    
signals:

    void params_changed(CParameters params);

private slots:

    void on_saveButton_clicked();

    void on_loadButton_clicked();

public slots:

    void on_applyButton_clicked();

private:
    Ui::Preferences *ui;

};

#endif // PREFERENCES_H
