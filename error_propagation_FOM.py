import math
frequency = (1./8.89244e-05) #[s-1]                                                 
inel = 8. * 1e-24




def peak_gauss(fit):
    peak = fit.GetMaximum()
    return peak

def propagate_errors_fom(calc, peak, func, params, errors,width_y,n1n2,fit):
    assert len(params) == len(errors)

    # Reset all function parameters to the
    # best fit values.
    for i in range(func.GetNpar()):
        func.SetParameter(i, params[i])

    calc_errs = []

    for i in range(func.GetNpar()):

            # Up variation: parameter + 1 sigma


        func.SetParameter(i, params[i] + errors[i])



        mu_up = calc(fit)



        peak_up = peak_gauss(fit)


        fom_up = peak_up /( frequency * n1n2* 1e22 / (mu_up * width_y  * 2.0 * math.pi)/ (frequency / inel ) )

           # Down variation: parameter - 1 sigma


        func.SetParameter(i, params[i] - errors[i])
        mu_down = calc(fit)
        peak_down = peak_gauss(fit)
        fom_down = peak_down /( frequency * n1n2* 1e22 / (mu_down * width_y  * 2.0 * math.pi)/ (frequency / inel ) )
        #           print fom_down,"fom down"
            # Reset parameter
        func.SetParameter(i, params[i])

            # Store symmetrised uncertainty


        calc_errs.append((fom_up - fom_down) / 2)

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


def calc_with_error_fom(calc,peak, func, fit_result,width_y,n1n2,fit):

    import numpy as np
    from math import sqrt


    errors =propagate_errors_fom(calc, peak, func, fit_result.GetParams(), fit_result.GetErrors(), width_y, n1n2,fit)


    errvec = np.array(errors)

    corr = correlation_matrix(func, fit_result)
    err2 = np.dot(errvec,np.dot(corr,errvec.T))

    return  err2


def propagate_errors_fom_y(calc, func, params, errors, width_x, peak_x, n1n2, fit):
    assert len(params) == len(errors)

    # Reset all function parameters to the
    # best fit values.
    for i in range(func.GetNpar()):
        func.SetParameter(i, params[i])

    calc_errs = []

    for i in range(func.GetNpar()):

            # Up variation: parameter + 1 sigma


        func.SetParameter(i, params[i] + errors[i])


        #print func.GetParameter([i]),"param"
        mu_up = calc(fit)
#        print mu_up,"mu up"
#        peak_up = peak(fit)
#            print peak_up,"peak up"
        fom_up = peak_x /( frequency * n1n2* 1e22 / (width_x * mu_up  * 2.0 * math.pi)/ (frequency / inel ) )
#            print fom_up,"fom up"
           # Down variation: parameter - 1 sigma


        func.SetParameter(i, params[i] - errors[i])
        mu_down = calc(fit)
#        peak_down = peak(fit)
        fom_down = peak_x /( frequency * n1n2* 1e22 / (width_x * mu_down  * 2.0 * math.pi)/ (frequency / inel ) )
        #           print fom_down,"fom down"
            # Reset parameter
        func.SetParameter(i, params[i])

            # Store symmetrised uncertainty


        calc_errs.append((fom_up - fom_down) / 2)

    return calc_errs





def calc_with_error_fom_y(calc, func, fit_result, width_x, peak_x, n1n2, fit):

    import numpy as np
    from math import sqrt


    errors =propagate_errors_fom_y(calc, func, fit_result.GetParams(), fit_result.GetErrors(), width_x,peak_x, n1n2,fit)


    errvec = np.array(errors)

    corr = correlation_matrix(func, fit_result)
    err2 = np.dot(errvec,np.dot(corr,errvec.T))

    return  err2

