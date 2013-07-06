#include "preferences.h"
#include "ui_preferences.h"

Preferences::Preferences(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Preferences) {

    ui->setupUi(this);

}

Preferences::~Preferences() {

    delete ui;

}

void Preferences::on_saveButton_clicked() {

}

void Preferences::on_loadButton_clicked() {

}

void Preferences::on_applyButton_clicked() {



}
