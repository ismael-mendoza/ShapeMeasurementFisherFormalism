import argparse
import defaults

def main():




    parser = argparse.ArgumentParser(description = 'Display different results and plots from the galaxies in a given file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('--wdir', default = 'output', metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('-f', '--filename', default = 'galaxies', metavar = 'FILENAME', type = str,
        help = 'Specify a file where galaxies will be registered or read from.')
    parser.add_argument('--plots-dir', default = 'plots', metavar = 'PLOTDIRECTORY', type = str,
    help = 'Specify a directory name where the plots will be saved')
    parser.add_argument('--hide', action = 'store_true',
    help = 'Do not show plots produced')
    parser.add_argument('--draw-galaxy', action = 'store_true',
    help = 'Show original galaxies')
    parser.add_argument('--partials', action = 'store_true',
        help = 'Show partial derivative images.')
    parser.add_argument('--second-partials', action = 'store_true',
        help = 'Show second partial derivative images.')
    parser.add_argument('--fisher', action = 'store_true',
        help = 'Show fisher matrix images.')
    parser.add_argument('--fisher-chi2', action = 'store_true',
        help = 'Show fisher matrix images with chi2 calculated values on top.')
    parser.add_argument('--covariance', action = 'store_true',
        help = 'Show covariance matrix elements')
    parser.add_argument('--correlation', action = 'store_true',
        help = 'Show correlation matrix elements.')   
    parser.add_argument('--bias-matrix', action = 'store_true',
        help = 'Show bias matrix images.')
    parser.add_argument('--biases', action = 'store_true',
        help = 'Show bias images.')
    parser.add_argument('--values', action = 'store_true', 
        help = 'Show values on top of appropiate images (fisher, biases).')


    #have to look and read the specified file. 
        #then based on the galaxy type and the parameters draw the galaxy.
        #each row could potentially represent a galaxy.




    #possibly add options to write some results down and compare with fits.py??? or useless 


plt.rc('text', usetex=False)

if args.hide:
    pass

if args.draw_galaxy:
    figure1, subplt= plt.subplots(1,1)
    figure1.suptitle('Initial Galaxy', fontsize = 20)
    drawPlot(subplt, gal_image.array)
    SaveFigureToPdfAndOpen(figure1, 'figure1.png')

if args.partials:
    figure2 = plt.figure() 
    figure2.suptitle('Derivatives of model with respect to each parameter', fontsize = 20)
    for i in range(num_params):
        ax = figure2.add_subplot(1,6,i+1)
        drawPlot(ax, Ds_gal[param_names[i]], title = param_names[i])
    SaveFigureToPdfAndOpen(figure2, 'figure2.png')




if args.fisher:
    figure3 = plt.figure()
    figure3.suptitle('Fisher matrix elements ', fontsize=14, fontweight='bold')
    for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure3.add_subplot(6,6, 6 * i + j + 1)
                    drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(param_names[i] )
                    if(i == len(Ds_gal) - 1):
                        ax.set_xlabel(param_names[j])
    SaveFigureToPdfAndOpen(figure3, 'figure3.png')


if args.fisher and args.value:
    figure4 = plt.figure()
    figure4.suptitle('Fisher matrix elements with values ', fontsize=14, fontweight='bold')
    for i in range(num_params): 
        for j in range(num_params):
            if(i >= j):
                ax = figure4.add_subplot(6,6, 6 * i + j + 1)
                drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                if(j == 0):
                    ax.set_ylabel(param_names[i])
                if(i == len(Ds_gal) - 1):
                    ax.set_xlabel(param_names[j])
                ax.text(-20,0, str(round(FisherM[param_names[i],param_names[j]],5)), fontweight='bold')
    SaveFigureToPdfAndOpen(figure4, 'figure4.png')

if args.fisher_chi2:
    figure5 = plt.figure()
    figure5.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=14, fontweight='bold')
    for i in range(num_params): 
        for j in range(num_params):
            if(i >= j):
                ax = figure5.add_subplot(6,6, 6 * i + j + 1)
                drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                if(j == 0):
                    ax.set_ylabel(param_names[i] )
                if(i == len(Ds_gal) - 1):
                    ax.set_xlabel(param_names[j])
                ax.text(-20,0, str(round(FisherM_chi2[param_names[i],param_names[j]],5)), fontweight='bold')
    SaveFigureToPdfAndOpen(figure5, 'figure5.png')

if args.covariance:
    figure6 = plt.figure()
    figure6.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
    for i in range(num_params): 
        for j in range(num_params):
            if(i >= j):
                ax = figure6.add_subplot(6,6, 6 * i + j + 1)
                drawPlot(ax, FisherM_images[param_names[i],param_names[j]] * 0)
                if(j == 0):
                    ax.set_ylabel(param_names[i])
                if(i == 6 - 1):
                    ax.set_xlabel(param_names[j])
                ax.text(-20,0, str(round(CovM[param_names[i],param_names[j]],5)), fontweight='bold')
    SaveFigureToPdfAndOpen(figure6, 'figure6.png')


########################change this to correlation and do the color thing to show when more correlated like David.
if args.correlation:
    figure7 = plt.figure()
    figure7.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
    for i in range(num_params): 
        for j in range(num_params):
            if(i >= j):
                ax = figure7.add_subplot(6,6, 6 * i + j + 1)
                drawPlot(ax, FisherM_images[param_names[i],param_names[j]] * 0)
                if(j == 0):
                    ax.set_ylabel(param_names[i])
                if(i == 6 - 1):
                    ax.set_xlabel(param_names[j])
                ax.text(-20,0, str(round(CovM[param_names[i],param_names[j]],5)), fontweight='bold') 
    SaveFigureToPdfAndOpen(figure7, 'figure7.png')

if args.second_partials:
    figure8 = plt.figure()
    figure8.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=14, fontweight='bold')
    for i in range(num_params): 
        for j in range(num_params):
            if(i >= j):
                ax = figure8.add_subplot(6,6, 6 * i + j + 1)
                drawPlot(ax, SecondDs_gal[param_names[i],param_names[j]])
                ax.text(0,20, "std:" + str(SecondDs_gal[param_names[i],param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
                if(j == 0):
                    ax.set_ylabel(param_names[i])
                if(i == len(Ds_gal) - 1):
                    ax.set_xlabel(param_names[j])

if args.bias_matrix:
    figuresOfBiasMatrix = [] 
    for i in range(num_params):
        figure = plt.figure()
        figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
        for j in range(num_params): 
            for k in range(num_params):
                if(j >= k):
                    ax = figure.add_subplot(6,6, 6 * j + k + 1)
                    drawPlot(ax, BiasM_images[param_names[i],param_names[j],param_names[k]]) 
                    if(k == 0):
                        ax.set_ylabel(param_names[j])
                    if(j == len(Ds_gal) - 1):
                        ax.set_xlabel(param_names[k])
        figuresOfBiasMatrix.append(figure)
    for i, figure in enumerate(figuresOfBiasMatrix): 
        SaveFigureToPdfAndOpen(figure, 'figure' + str(8) + '_' + str(i) + '.png') 

if args.bias_matrix and args.values: 
    figuresOfBiasMatrixNumbers = []    
    for i in range(num_params):
        figure = plt.figure()
        figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
        for j in range(num_params): 
            for k in range(num_params):
                if(j >= k):
                    ax = figure.add_subplot(6,6, 6 * j + k + 1)
                    drawPlot(ax, BiasM_images[param_names[i],param_names[j],param_names[k]]) 
                    if(k == 0):
                        ax.set_ylabel(param_names[j] )
                    if(j == len(Ds_gal) - 1):
                        ax.set_xlabel(param_names[k])
                    ax.text(-20,0, str(round(BiasM[param_names[i],param_names[j],param_names[k]],5)), fontweight='bold')
        figuresOfBiasMatrixNumbers.append(figure)

    for i, figure in enumerate(figuresOfBiasMatrixNumbers): 
        SaveFigureToPdfAndOpen(figure, 'figure' + str(9) + '_' + str(i) + '.png') 


if args.biases:
    figure10 = plt.figure() 
    figure10.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 15)
    for i in range(num_params):
        ax = figure10.add_subplot(1,6,i+1)
        drawPlot(ax, bias_images[param_names[i]], title = param_names[i])

    SaveFigureToPdfAndOpen(figure10, 'figure10.png')
 
 if args.biases and args.values:
    figure11 = plt.figure() 
    figure11.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 14)
    for i in range(num_params):
        ax = figure11.add_subplot(1,6,i+1)
        drawPlot(ax, bias_images[param_names[i]], title = param_names[i])
        ax.text(-20,0, str(round(biases[param_names[i]],5)), fontweight='bold')

    SaveFigureToPdfAndOpen(figure11, 'figure11.png')





    list_of_plots = np.array([int(n) for n in argv[2:]])

    #here we slice arg[2:] because we have to have a numpy array only with numbers for this to work and only the arguments after 'plot' are numbers.
    if((list_of_plots == 1).any()): 
 #this will open for preview because it is the default defined in my mac


    if((list_of_plots == 2).any()):  
        


    if((list_of_plots == 3).any()):



    if((list_of_plots == 4).any()):



    if((list_of_plots == 5).any()):


        

    if((list_of_plots == 6).any()):



    if((list_of_plots == 7).any()):


        

    if((list_of_plots == 8).any()):




    if((list_of_plots == 9).any()):

    if((list_of_plots == 10).any()):


    if((list_of_plots == 11).any()):



if __name__ == '__main__':
    main()