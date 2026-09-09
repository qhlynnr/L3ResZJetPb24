#include <iostream>
#include <string>
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TH1.h"

void CheckHists(TDirectory *dir) {
    TIter next(dir->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();

        if (obj->InheritsFrom("TDirectory")) {
            CheckHists((TDirectory*)obj);
        }
        else if (obj->InheritsFrom("TH1")) {
            TH1 *h = (TH1*)obj;

            cout << dir->GetPath() << "/" << h->GetName()
                 << " : entries = " << h->GetEntries()
                 << ", integral = " << h->Integral()
                 << endl;
        }
    }
}

void CheckEmptyHist(string filename) {
    TFile *f = TFile::Open(filename.c_str());

    if (!f || f->IsZombie()) {
        cout << "Could not open " << filename << endl;
        return;
    }

    CheckHists(f);

    f->Close();
}
