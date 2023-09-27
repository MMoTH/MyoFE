import numpy as np
import scipy.interpolate as interpol


def move_cb_distributions(self, delta_hsl):
    """ Moves cb distributions """

    delta_x = delta_hsl * self.implementation['filament_compliance_factor']

    if (self.implementation['kinetic_scheme'] == '3_state_with_SRX') or \
        (self.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):
        interp_positions = self.x - delta_x
        bin_indices = 2 + np.arange(0, self.no_of_x_bins)
        
        # Count bridges before
        before_heads = np.sum(self.y[bin_indices])
        self.y[bin_indices] = interpol.interp1d(self.x,
                                  self.y[bin_indices],
                                  kind='quadratic',
                                  fill_value=0,
                                  bounds_error=False)(interp_positions)
        #print 'y_vec in move cb'
        #print self.y

        # Make sure we don't have any negative populations
        temp = self.y[bin_indices]
        temp[np.nonzero(temp < 0.0)] = 0.0
        self.y[bin_indices] = temp

        # Count bridges now
        after_heads = np.sum(self.y[bin_indices])
        # These appear in M_on
        self.y[1] = self.y[1] + (before_heads - after_heads)

    if (self.implementation['kinetic_scheme'] == '4_state_with_SRX') or \
        (self.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):
        interp_positions = self.x - delta_x

        for i in range(0, 2):
            if (i == 0):
                ind = 2 + np.arange(0, self.no_of_x_bins)
            else:
                ind = 2 + self.no_of_x_bins + np.arange(0, self.no_of_x_bins)

            # Count bridges before
            before_heads = np.sum(self.y[ind])

            # Interpolation
            self.y[ind] = interpol.interp1d(
                            self.x,
                            self.y[ind],
                            kind='quadratic',
                            fill_value=0,
                            bounds_error=False)(interp_positions)

            # Make sure we don't have any negative populations
            temp = self.y[ind]
            temp[np.nonzero(temp < 0.0)] = 0.0
            self.y[ind] = temp

            # Count bridges now
            after_heads = np.sum(self.y[ind])
            # These appear in M_on
            self.y[1] = self.y[1] + (before_heads - after_heads)
