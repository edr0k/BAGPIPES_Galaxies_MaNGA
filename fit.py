import numpy as np
import pandas as pd
import bagpipes as pipes
from astropy.table import Table, join
import math

# Essa é a biblioteca do BAGPIPES, para baixar basta acessar: https://bagpipes.readthedocs.io/en/latest/index.html e seguir os
# os passos. O PyMultiNest é essencial e também precisa ser instalado, caso contrário o código não vai rodar, mas nesse
# mesmo link ele indica como fazer. Qualquer coisa pode me mandar email que eu te ajudo.


# A ordem dos filtros importa para o código, então tem que colocar na mesma ordem que a fotometria foi passada.
filter_list_splus = ['filters/uJAVA.dat',
                     'filters/F0378.dat',
                     'filters/F0395.dat',
                     'filters/F0410.dat',
                     'filters/F0430.dat',
                     'filters/gSDSS.dat',
                     'filters/F0515.dat',
                     'filters/rSDSS.dat',
                     'filters/F0660.dat',
                     'filters/iSDSS.dat',
                     'filters/F0861.dat',
                     'filters/zSDSS.dat']


mags = ["u_petro", "J0378_petro", "J0395_petro", "J0410_petro", "J0430_petro", "g_petro", "J0515_petro", "r_petro",
        "J0660_petro", "i_petro", "J0861_petro", "z_petro"]
mags_err = ["e_u_petro", "e_J0378_petro", "e_J0395_petro", "e_J0410_petro", "e_J0430_petro", "e_g_petro",
            "e_J0515_petro","e_r_petro", "e_J0660_petro","e_i_petro", "e_J0861_petro", "e_z_petro"]
#####################################################################################################################
# Dicionários de parâmetros fixos e livres, nessa parte existe uma boa flexibilidade do fitting.
# Segue o link da documentação do BAGPIPES, com alguns exemplos de como usar:
# https://bagpipes.readthedocs.io/en/latest/fit_instructions.html

# A primeiro momento estou usando um modelo simples de Burst, deixando a idade, metalicidade, massa de formação para
# serem ajustados. Nunca tinha usado ele e descobri que se eu faço só isso, ele não faz muito bem a SFH e não sei
# muito bem o porquê.
# defining fit parameteres
delayed = {}
delayed["age"] = (0., 15.)
delayed["tau"] = (0., 15.)
delayed["massformed"] = (1., 15.)  # 10
delayed["metallicity"] = (0., 2.5)  # 0.02

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0., 2.)
dust["eta"] = 2.

nebular = {}
nebular["logU"] = -3.

fit_instructions = {}
#fit_instructions["redshift"] = (0.,0.05)
fit_instructions["exponential"] = delayed
fit_instructions["dust"] = dust
fit_instructions["nebular"] = nebular
#fit_instructions["veldisp"] = (0.,400.)
df = pd.read_csv('galaxies_idr4.csv')

#Função load data necessária para o fit. Ela pega a fotometria do ID que é passado e transforma para fluxo em mJy,
# assim como a incerteza. A conta do fluxo fui eu quem fez, espero estar certa. Se não estiver, por favor me avise!
def load_data_splus(id):
        print(id)
        # try:
        #     object, field = id.split('_')
        # except:
        #     object1, object2, field = id.split('_')
        #     object = object1 + '_' + object2

        fluxes = []
        fluxerrs = []
        galaxy_param = df[(df['ID'] == id)]
        #print(galaxy_param[mags])
        #print(galaxy_param[mags_err])

        for k in range(0, len(mags)):
            m = galaxy_param[mags[k]]
            if (math.isnan(m.values)) | (m.values == np.inf) | (m.values == 99.) | (m.values == -99.):
                f = np.array([99.])
                delta_f = np.array([99.])
            else:
                f = 10**(9.56) * 10**(-m/2.5)  # flux in mJy
                delta_m = galaxy_param[mags_err[k]]
                delta_f = f * (1/2.5) * np.log(10) * delta_m

            #print(m, delta_m, f, delta_f)

            fluxes.append(f)  #mJy
            fluxerrs.append(delta_f)  #mJy
        photometry = np.c_[fluxes, fluxerrs]

        print(photometry)
        return photometry

IDs = df['ID']
redshifts = df['zml']
print(redshifts)
fit_cat = pipes.fit_catalogue(IDs, fit_instructions, load_data_splus, spectrum_exists=False,
                              cat_filt_list=filter_list_splus, run="galaxies_idr4", make_plots=True,
                              full_catalogue=True, redshifts=redshifts)

fit_cat.fit(verbose=False)
