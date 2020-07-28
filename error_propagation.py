

import math
frequency = (1./8.89244e-05) #[s-1]  
inel = 8. * 1e-24




def width(fit):

        # maximum of the fit fucntion                                                                                              
    peak = fit.GetMaximum()
#        print peak,"peak compare me"                                                                                              
        # position of the peak                                                                                                     


        # limits for the integration                                                                                               

    intlimit = 0.1

        # integral of the fit function( it gives me area under the curve)                                                          

    integral = fit.Integral(-intlimit,intlimit)

        # formula for finding the cap sigma                                                                                        
    width_beam =(1 / math.sqrt (2 * math.pi)) * integral / peak

    return width_beam



def peak_gauss(fit):
    peak = fit.GetMaximum()
    return peak














def propagate_errors(calc, func, params, errors,fit):

    assert len(params) == len(errors)

    # Reset all function parameters to the
    # best fit values.
    for i in range(func.GetNpar()):
        func.SetParameter(i, params[i])

    calc_errs = []

    for i in range(func.GetNpar()):
        # Up variation: parameter + 1 sigma
        
        func.SetParameter(i, params[i] + errors[i])

        res_up = calc(fit)

#        print res_up,"up"
        # Down variation: parameter - 1 sigma
        negative=func.SetParameter(i, params[i] - errors[i])


        res_down = calc(fit)
#        print res_down,"down"

        # Reset parameter
        func.SetParameter(i, params[i])

        # Store symmetrised uncertainty
        residual=calc_errs.append((res_up - res_down) / 2)
  #  print calc_errs,"residual"


    return calc_errs




def correlation_matrix(func, fit_result):
    import numpy as np

    n = func.GetNpar()
    corrmatrix = fit_result.GetCorrelationMatrix()
    corr = np.ndarray((n,n))
    for i in xrange(n):
        for j in xrange(n):
            corr[i,j] = corrmatrix(i,j)

    return corr






def propagate_errors_lumi(calc, func, params, errors,width_y,n1n2,fit):

    assert len(params) == len(errors)

    # Reset all function parameters to the
    # best fit values.
    for i in range(func.GetNpar()):
        func.SetParameter(i, params[i])

    calc_errs = []

    for i in range(func.GetNpar()):
        # Up variation: parameter + 1 sigma
        
        func.SetParameter(i, params[i] + errors[i])

        res_up = calc(fit)
        mu_up = frequency * n1n2* 1e22 / (res_up * width_y  * 2.0 * math.pi)/ (frequency / inel )
        




#        print res_up,"up"
        # Down variation: parameter - 1 sigma
        func.SetParameter(i, params[i] - errors[i])


        res_down = calc(fit)
        mu_down = frequency * n1n2 * 1e22 / (res_down * width_y  * 2.0 * math.pi)/ (frequency / inel )


#        print res_down,"down"

        # Reset parameter
        func.SetParameter(i, params[i])

        # Store symmetrised uncertainty
        residual=calc_errs.append((mu_up - mu_down) / 2)
  #  print calc_errs,"residual"


    return calc_errs






def calc_with_error(calc, func, fit_result,fit):
    import numpy as np
    from math import sqrt

#    nom = calc()

    errors = propagate_errors(calc, func, fit_result.GetParams(), fit_result.GetErrors(),fit)
    errvec = np.array(errors)
    #print errvec,"errvec"
    corr = correlation_matrix(func, fit_result)
    err2 = np.dot(errvec,np.dot(corr,errvec.T))

#    print sqrt(err2),"err2"
    
    return sqrt(err2)









def calc_with_error_lumi(calc, func, fit_result,width_y,n1n2,fit):
    import numpy as np
    from math import sqrt

#    nom = calc()

    errors = propagate_errors_lumi(calc, func, fit_result.GetParams(), fit_result.GetErrors(),width_y,n1n2,fit)
    errvec = np.array(errors)
    #print errvec,"errvec"
    corr = correlation_matrix(func, fit_result)
    err2 = np.dot(errvec,np.dot(corr,errvec.T))

#    print sqrt(err2),"err2"
    
    return err2


