import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict
from sklearn.metrics import matthews_corrcoef 
from matplotlib.colors import ListedColormap
from plot_graph import test
from matplotlib.lines import Line2D
import scipy.stats as st
from sklearn.svm import SVC
from bioservices.kegg import KEGG
import numpy as np
from matplotlib import rc
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from met_names import get_reaction_names
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx
from igraph import *
from colormap import rgb2hex
from scipy.stats import linregress
rc('font',**{'family':'serif','serif':['STIX']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 36})
plt.rcParams['figure.figsize'] = 12,10


def process():
	preds = pd.read_pickle('all_preds.pickle')
	#preds = preds.drop(columns = ['Apts_Ser', 'LIPASE_METHYL_OCDCEA_MG327'])
	accs = pd.read_pickle('all_accs_second.pickle')
	accm = accs.max(axis = 1)
	acc_bad = accm[accm < 0.7].index
	preds = preds.drop(columns = acc_bad)
	#preds = preds.sort_values(by = 'ind')
	preds_inds = preds['ind']
	#preds = preds.drop(columns = ['ind'])
	#d = {0:'Septum', 1:'Protein', 2:'Metabolism', 3:'DNA', 4:'Non-essential', 5:'RNA', 6:'Slow-growing', 7:'Known Error'}
	d = {}
	labels = pd.read_pickle('labels.pickle')
	labels = labels.dropna()
	full = labels['Class']
	pca = PCA(n_components = 2)
	out = pd.DataFrame(np.hstack((np.expand_dims(preds_inds, axis = 1), pca.fit_transform(preds[preds.columns[:-1]]))), columns = ['reps', 'Principle Component 1', 'Principle Component 2'])
	l = [0] * len(full)
	j = 0
	sns_hue = []
	for i in set(full):
		d[j] = i
		j = j+1
		col = np.random.choice(range(256), size = 3)/256
		sns_hue.append(col)
		inds = np.where(full == i)
		for x,y in zip(inds[0], [col]*len(inds[0])):
			l[x] = y
	labels['colours'] = l
	new = labels.merge(out, on = 'reps')
	preds_keep = np.arange(preds_inds.shape[0])[np.in1d(preds_inds.values, new['reps'].values)]
	preds = preds.loc[preds_keep]
	new = new.drop_duplicates(['reps'])
	new = new.sort_values(by = 'reps')
	preds = preds.sort_values(by = 'ind')
	preds = preds.drop(columns = 'ind')
	preds = preds.reset_index()
	preds = preds.drop(columns = 'index')
	new = new.reset_index()
	new = new.drop(columns = 'index')
	return new, preds


def multiassign(d, keys, value):
        for k in keys:
                d[k] = value
        return d

def plot_regression():
	fluxs = pd.read_pickle('/home/sl13479/Documents/PhD_Docs/Oli_Suite/Testing/output/wild_test/1/fluxs_wild_test_1.pickle')
	fluxs[3].plot()
	reg = linregress(fluxs.index.values, fluxs[3].values)
	plt.plot(fluxs.index.values, ((fluxs.index.values) * reg.slope) + reg.intercept)
	plt.xlabel('Time (s)')
	plt.ylabel('Flux')
	plt.savefig('/figs/regression.pdf')
	return

def plot(new):
	#new.plot.scatter('Principle Component 1', 'Principle Component 2', c = new['colours'], figsize =[15,10])
	col = []
	for idx, val in enumerate(set(new['Class'])):
		col.append(np.random.choice(range(256), size = 3)/256)
	ax = sns.scatterplot(new['Principle Component 1'], new['Principle Component 2'], hue = new['Class'], palette = col, s = 100)
	legend_dict = {}
	patch_list = []
	#for n in set(full):
	#	name = mpatches.Patch(color = df['colours'][np.where(full == n)[0][0]], label = n)
	#	legend_dict['patch_' + str(n)] = name
	#	patch_list.append(name)
		#plt.legend(handles=patch_list)
	plt.savefig('./figs/pca_labelled1.svg')
	plt.savefig('./figs/pca_labelled1.pdf')
	return ax

def spectral():
	preds = pd.read_pickle('all_preds_full.pickle')
	pca = PCA(n_components=2)
	out = pca.fit_transform(preds)
	#inds = np.where(preds['Apts_Ser'] == 1)
	col_one = np.random.choice(range(256), size = 3)/256
	l = [col_one] * len(preds)
	col_two = np.random.choice(range(256), size = 3)/256
	for x, y in zip(inds[0], [col_two]*len(inds[0])):
		l[x] = y
	pca_df = pd.DataFrame(out, columns = ['Principle Component 1', 'Principle Component 2'])
	pca_df['colours'] = l
	pca1 = pd.DataFrame([pca.components_[0], pca.components_[1]], index = preds.columns)
	#pca_df.plot.scatter('Principle Component 1', 'Principle Component 2', c = pca_df['colours'], figsize =[15,10])
	#plt.savefig('/figs/pca_apt_ser.pdf')
	#plt.clf()
	#plt.plot(pca.components_[0], 'ko')
	#plt.plot(pca.components_[1], 'co')
	#plt.savefig('/figs/eigs.pdf')
	#plt.clf()
	#sns.heatmap(pca.get_covariance())
	#plt.savefig('/figs/covariance.pdf')
	return

def drivers():
	stoich = np.genfromtxt('full_stoich.csv', delimiter = ',')
	reac_adj = np.matmul(np.transpose(stoich), stoich)
	reac_graph = nx.from_numpy_matrix(reac_adj)
	matching = nx.maximal_matching(reac_graph)
	nodes = list(reac_graph.nodes)
	driver_nodes = list(set(nodes).symmetric_difference([item for t in matching for item in t]))
	return driver_nodes, reac_adj

def graph(drivers):
	stoich = np.genfromtxt('full_stoich.csv', delimiter = ',')
	hubs = np.where(stoich.sum(axis = 1) > 10)
	#stoich[hubs,:] = 0
	reac_adj = np.matmul(np.transpose(stoich), stoich)
	reac_graph = nx.from_numpy_matrix(reac_adj)
	G = Graph()
	G.add_vertices(reac_graph.nodes)
	G.add_edges(reac_graph.edges)
	col = [0] * len(G.vs.indices)
	for n in G.vs.indices:
		col[n] = 'rgb' + str('(100, 100, 100)')
	for x in drivers:
		col[x] = 'rgb' + str('(40, 10, 10)')
	G.vs['color'] = col
	return G

def define_drivers(driver_nodes, reac_adj):
	driver_dict = {}
	for x in driver_nodes:
		neighbours = np.where(reac_adj[:, x] != 0)
		driver_dict[x] = neighbours[0]
	return driver_dict

def get_pathways(EC, k, myco_names):
	gene = k.find('mge', 'ec:' + str(EC))
	pathways = None
	name = None
	if len(gene) == 1:
		for x in myco_names:
			gene = k.find(x, 'ec:' + str(EC))
			if len(gene) > 1:
				print(x)
				name = gene.split('\t')[0].split(':')[1]
				pathways = k.get_pathway_by_gene(name, x)
				break
	elif len(gene) > 1:
		x = 'mge'
		name = gene.split('\t')[0].split(':')[1]
		pathways = k.get_pathway_by_gene(name, 'mge')
		#print(driver_dict[key])
		print(pathways)
	else:
		pathways = None
	return pathways, x, name

def get_driver_pathways(driver_dict):
	l_path = []
	l_org = []
	l_gene = []
	d_node = []
	node = []
	mapping = pd.read_csv('/home/sl13479/Documents/PhD_Docs/metabolism_mapping/reaction_gene_mapping.csv')
	k = KEGG()
	orgs = k.lookfor_organism('mycoplasma')
	myco_names = [x.split(' ')[1] for x in orgs]
	for key in driver_dict:
		main_enzyme = mapping['EC'][[key]]
		if isinstance(main_enzyme, str):
			enzymes = main_enzyme
			node.append(key)
			#d_node.append(key)
		else:
			enzymes = mapping.iloc[driver_dict[key]]['EC'].dropna().drop_duplicates()
			node.append(enzymes.index.values)
		for x in enzymes:
			d_node.append(key)
			pathways, org, gene  = get_pathways(x, k, myco_names)
			l_path.append(pathways)
			l_org.append(org)
			l_gene.append(gene)
			print(x)
	info = pd.DataFrame(np.column_stack([l_path, l_org, l_gene, d_node, np.concatenate(node).ravel()]), columns = ['pathways', 'organism', 'gene', 'driver_node', 'node'])
	info = info.dropna()
	return info

def class_centroid(new):
	all_means = []
	classes = []
	for n in set(new['Class']):
		mean = new[new['Class'] == n][['Principle Component 1', 'Principle Component 2']].mean().values 
		all_means.append(mean)
		classes.append(n)
	means_df = pd.DataFrame(all_means, index = classes)
	means_df = means_df.drop('Known Error')
	return means_df

def svm(df, ax, classif, means_df, name):
	lim = int(np.rint(0.8*len(df)))
	train_df = df.iloc[0:lim]
	test_df = df.iloc[lim+1:len(df)]
	train_class = classif[0:lim]
	test_class = classif[lim+1:]
	clf = SVC(kernel = 'linear')
	clf.fit(train_df, train_class)
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	acc = clf.score(test_df, test_class)
	pred = clf.predict(np.concatenate(means_df.values).reshape([7,2]))
	x = np.linspace(xlim[0], xlim[1] , 30)
	y = np.linspace(ylim[0], ylim[1], 30)
	Y, X = np.meshgrid(y, x)
	xy = np.vstack([X.ravel(), Y.ravel()]).T
	P = clf.decision_function(xy).reshape(X.shape)
	ax.contour(X, Y, P, colors='k', levels=[0], alpha=0.5, linestyles=['-'])
	#x = [d[0] for d in means_df['means'].values]
	#y = [d[1] for d in means_df['means'].values]
	#plt.plot(x,y, 'ko')
	plt.savefig('./figs/pca/pca_' + str(name) + '.svg', bbox_inches = 'tight')
	plt.savefig('./figs/pca/pca_' + str(name) + '.pdf', bbox_inches = 'tight')
	plt.savefig('./figs/pca/pca_' + str(name) + '.eps', bbox_inches = 'tight')
	return acc, pred

def get_reacs_inds(sim, preds, flux_names):
	if isinstance(sim, int):
		x = preds.iloc[sim]
	else:
		x = sim
	col_names = x[x != 0].index.values
	true_inds = [flux_names.index(x) for x in col_names]
	col_ind = [x.index.get_loc(y) for y in x[x != 0].index.values]
	return true_inds, col_ind

def groups(preds, df, clust, thresh):
	inds = df[df['Class'] == clust].index.values
	weights = preds.iloc[inds].sum()
	weights[weights < ((thresh/80) * len(inds))] = 0 #find noise threshold from sd of wildtype
	norm = weights/len(inds)
	col = np.random.choice(range(128), size = 3)
	full_col = norm.apply(lambda x: col + (np.array([255,255,255]) - col)*(1-x))
	#col_list = [tuple(x) for x in full_col.values]
	per = max(norm)
	return full_col, norm, per

def colourbar(full_col):
	hexc = np.unique([rgb2hex(int(x[0]), int(x[1]), int(x[2])) for x in full_col.values])
	cmap = ListedColormap(sns.color_palette(np.unique(hexc)), name = 'cust')
	sm = plt.cm.ScalarMappable(cmap=cmap.reversed())
	#cbar = plt.colorbar(sm, orientation='horizontal')
	return sm, cmap
	

def fit_expon():
	preds = pd.read_pickle('preds_wild.pickle')
	smean = preds.sum().values.mean()
	rate = 1./smean
	x = np.linspace(0, preds.sum().max(), preds.sum().max())
	dist_exp = st.expon.pdf(x, scale = 1./rate)
	point = st.expon.interval(0.95, loc=0, scale = 1./rate)
	fig, ax = plt.subplots()
	sns.distplot(preds.sum().values, kde = False, ax=ax)
	ax.plot(x, dist_exp)
	return point[1]

def plot_graph(inds, colour_grade):
        stoich = np.genfromtxt('full_stoich.csv', delimiter = ',')
        hubs = np.where(stoich.sum(axis = 1) > 10)
        stoich[hubs,:] = 0
        reac_adj = np.matmul(np.transpose(stoich), stoich)    
        G = nx.from_numpy_matrix(reac_adj)
        reac_g = Graph()
        reac_g.add_vertices(G.nodes)
        reac_g.add_edges(G.edges)
        d = {}
        col = multiassign(d, [x for x in range(645)], [255, 255, 255])
        if colour_grade is None:
                col = multiassign(col, [x for x in inds], 'cadetblue')
                #col[85] = 'crimson'
        else:
                for idx, val in enumerate(inds):
                        col[val] = 'rgb' + str(colour_grade[idx])
        reac_g.vs['color'] = [col[x] for x in reac_g.vs.indices]
        ids = reac_g.vs.select(color_eq=[255, 255, 255])
        g = reac_g.copy()
        g.delete_vertices(reac_g.vs(color_ne = [255, 255, 255]))
        reac_g.delete_vertices(ids)
        #layout = reac_g.layout_kamada_kawai()
        #plot(reac_g, layout = layout)
        return reac_g, col, g, ids

		
def plot_reacs(new, preds, name):
	df = new
	fig, ax = plt.subplots()
	col_one = np.random.choice(range(128), size = 3)/256
	l = [col_one] * len(df)
	col_two = np.random.choice(range(128, 256), size = 3)/256
	inds = np.where(preds[name] == 1)
	for x,y in zip(inds[0], [col_two]*len(inds[0])):
		l[x] = y
		#l[x] = sum(y)
	df['colours'] = l
	colours = [col_one, col_two]
	g = sns.scatterplot(new['Principle Component 1'], new['Principle Component 2'], style = preds[name].values, markers = ['o', 'D'], hue = preds[name].values, palette = colours[0:len(np.unique(preds[name].values))], s = 100, linewidth=0.7, ax=ax) #hue = preds[name].values, s = 40)
	custom = [Line2D([], [], marker='o', color=col_one, linestyle='None'),
		  Line2D([], [], marker='D', color=col_two, linestyle='None')]
	plt.legend(custom, ['Normal', 'Abnormal'])
	plt.savefig('./figs/pca/pca_' + str(name) + '.pdf')
	plt.savefig('./figs/pca/pca_' + str(name) + '.svg')
	plt.savefig('./figs/pca/pca_' + str(name) + '.eps')
#sum_l = np.array([sum(x) for x in l])
	#c = np.zeros(len(l))
	#for idx, val in enumerate(np.unique(sum_l)):
	#	c[np.where(sum_l == val)] = idx
	c = preds[name].values
	return ax, c

def mutual_info(new, preds, classif):
	boo = new['Class'] == classif
	labels = boo.apply(lambda x: int(x))
	I = preds.apply(metrics.adjusted_mutual_info_score, args=[labels])
	sns.scatterplot(x = np.linspace(0, len(I), len(I)), y = I.values)
	plt.show()
	plt.clf()
	return I

def entropy(new, preds):
	probs = np.zeros((len(preds.columns), len(np.unique(new['Class']))))
	i = 0
	for n in np.unique(new['Class']):
		inds = new[new['Class'] == n].index.values
		norm = (preds.iloc[inds].sum())/len(inds)
		probs[:,i] = norm.values
		i = i+1
	p = pd.DataFrame(probs, columns = np.unique(new['Class']))
	H = p.apply(lambda x: -1 * x*np.log(x)) #.sum(axis = 1)
	#H = H[H < -0.2]
	
	return probs, H

def driver_modes(new, preds):
	driver = pd.read_pickle('driver_df.pickle')
	names = get_reaction_names()
	#driver_names = []
	#for x in np.unique(driver['driver_node']):
	#	if names[x] in preds.columns:
	#		driver_names.append(names[x])	
	#preds_d = preds[driver_names]
	driver_name = []
	corr = []
	#corr_d = []
	for d in np.unique(driver['driver_node']):
		#corr = []
		if names[d] in preds.columns:
			corr.append([matthews_corrcoef(preds[x], preds[names[d]]) for x in preds.columns])		
	#		corr_d.append([matthews_corrcoef(preds_d[x], preds_d[names[d]]) for x in preds_d.columns])
	
			driver_name.append(names[d])
	out = pd.DataFrame(np.transpose(corr), columns = driver_name) 
	#out_d = pd.DataFrame(corr_d, columns = driver_names)
	clustering = AgglomerativeClustering(n_clusters = 6).fit(preds[driver_names].transpose())
	clust_dict = defaultdict(list)
	for idx, val in enumerate(clustering.labels_):
		clust_dict[val].append(driver_name[idx])
	return out, clust_dict

def sort(out, clust_dict):
	less = out[out > 0.9].dropna(how = 'all')
	d = defaultdict(list)
	m = less.max(axis = 1)
	for idx, row in less.iterrows():
		name = out.columns[np.where(m.loc[idx] == row)].values
		if len(name) == 1:
			d[name[0]].append(idx)
	new_d = defaultdict(list)
	for key in clust_dict:
		for val in clust_dict[key]:
			[new_d[key].append(x) for x in d[val]]
	return d, new_d

def drive_plot(new, preds, new_d):
	for key in new_d:
		if len(np.unique(new_d[key])) > 1:
			reas = preds.iloc[:, new_d[key]].sum(axis = 1)
			rea_norm = reas/reas.max()
			col = np.random.choice(range(199), size = 3)
			full_col = rea_norm.apply(lambda x: col + (np.array([200,200,200]) - col)*(1-x))
		#hexc = [rgb2hex(int(x[0]), int(x[1]), int(x[2])) for x in full_col.values]  
		#cmap = ListedColormap(sns.color_palette(np.unique(hexc)), name = 'cust')
			sm, cmap = colourbar(full_col)
		#ax = plt.scatter(x = new['Principle Component 1'], y = new['Principle Component 2'], c = hexc)
			fig, ax = plt.subplots()
			sns.scatterplot(x = new['Principle Component 1'], y = new['Principle Component 2'], palette = cmap, hue = rea_norm, s = 80, legend = False, ax=ax, linewidth = 0.4)
			cbar = plt.colorbar(sm)
			plt.savefig('./figs/drivers/driver_' + str(key) + '.pdf')
			plt.savefig('./figs/drivers/driver_' + str(key) + '.svg')
			plt.savefig('./figs/drivers/driver_' + str(key) + '.eps')
	return

def single(new, preds, name):
	c = new[new['Class'] == name]
	c_preds = preds.iloc[c.index.values]
	pca = PCA(n_components = 2)
	out = pca.fit_transform(c_preds)
	c = c.drop(columns = ['Principle Component 1', 'Principle Component 2'])
	c['Principle Component 1'] = out[:,0]
	c['Principle Component 2'] = out[:,1]
	cov = pca.get_covariance()
	return c, c_preds, cov

def group(new, preds, name, H):
	cut = H[name].dropna().mean()
	inds = H[name][H[name] < cut].index.values
	corr = []
	c = [int(x) for x in new['Class'] == name]  
	for col in preds.columns:
		corr.append(matthews_corrcoef(c, preds[col]))
	return pd.Series(corr)

def intersect(inds, inds_wt):
	#e.g. inds = [1343, 2282, 1582]
	#     inds = [2088, 2738, 2876, 3278]
	wt = set(np.where(preds.iloc[inds_wt].sum() > 0)[0])
	a = set(np.where(preds.iloc[inds[0]])[0]).difference(wt)
	b = set(np.where(preds.iloc[inds[1]])[0]).difference(wt)
	c = set(np.where(preds.iloc[inds[2]])[0]).difference(wt)
	abc = (a.intersection(b)).intersection(c)
	ab = (a.intersection(b)).symmetric_difference(abc)
	ac = (a.intersection(c)).symmetric_difference(abc)
	a_ = a.difference(b.union(c))
	b_ = b.difference(a.union(c))
	c_ = c.difference(b.union(a))
	return a_, ac, ab, abc

def find_points(new, means_df):
	dist = pd.DataFrame(cdist(new[['Principle Component 1', 'Principle Component 2']], np.expand_dims(means_df.mean(), axis = 0))).sort_values(by=0).iloc[0:10]
	inds = dist.index.values
	m = new.iloc[inds][['Principle Component 1', 'Principle Component 2']]
	d = cdist(m, m)
	first = np.argmin(d[np.nonzero(d)])
	return


#def dist(df):
#	d = cdist(df, df)
#	d_mean = df - df.mean()
#	d_euclid = d_mean.pow(2).sum(axis = 1).apply(np.sqrt)
#	first = np.where(d_euclid == d_euclid.max())
#	second = np.argmax(d[first])
#	d[first] = d[second] + d[first]
#	d = np.delete(d, second, axis = 0)
#	d = np.delete(d, second, axis = 1)
#	third = np.argmax(d[first])
#	inds = df.iloc[[first[0][0], second, third].index.values
#	return inds

def view(n, new, means_df):	
	df = new[new['Class'] == n][['Principle Component 1', 'Principle Component 2']]
	d_mean = df - df.mean()
	d_euclid = d_mean.pow(2).sum(axis = 1).apply(np.sqrt)
	first = df.iloc[np.where(d_euclid == d_euclid.max())]
	means_df = means_df.drop(index = n)
	other = np.argmin(cdist(means_df, first))
	#name = means_df.iloc[other].name
	close = np.argmin(cdist(new[new['Class'] != n][['Principle Component 1', 'Principle Component 2']], first))
	last = new[new['Class'] != n].iloc[close]
	wild = np.argmin(cdist(new[new['Class'] == 'Non Essential'][['Principle Component 1', 'Principle Component 2']], first))
	wild_ind = new[new['Class'] == 'Non Essential'].iloc[wild]
	return last, new.iloc[first.index.values[0]], wild_ind

	
	
if __name__ == '__main__':
	new, preds = process()
	#get processed data (neural network predictions for all simulations, and labels/pca)
	flux_names = get_reaction_names()
	threshold = fit_expon()
	col, norm, per = groups(preds, new, 'Protein', threshold)
	inds, col_ind = get_reacs_inds(norm, preds, flux_names)
	g, cols, G, ids = plot_graph(inds, [tuple(x) for x in col.iloc[col_ind].values])
	sm, cmap = colourbar(col.values)
	test(g, sm, mut)
	layout = g.layout_kamada_kawai()
	plot(g, layout = layout)
	for mut in np.unique(new['Class']):
		col, norm, per = pf.groups(preds, new, mut, threshold)
		inds, col_ind = pf.get_reacs_inds(norm, preds, flux_names)
		g, cols, G, ids = pf.plot_graph(inds, [tuple(x) for x in col.iloc[col_ind].values])
		sm, cmap = pf.colourbar(col)
		pg.test(g, sm, mut) #(import plot_full as pf, import plot_graph as pg)

	names = flux_names
#spectral()
#plot()
	all_accs = []
	means_df = class_centroid(new) 
	driver_nodes, reac_adj = drivers()
	driver_names = [names[x] for x in driver_nodes]
	driver_dict = define_drivers(driver_nodes, reac_adj)
	info = get_driver_pathways(driver_dict)
	driver_df = pd.read_pickle('driver_df.pickle')
	all_pred = []
	for n in info['node']:
		if names[n] in preds.columns.values:
			ax, c = plot_reacs(new, preds, names[n])
			if len(np.unique(c)) > 1:
				acc, pred = svm(new[['Principle Component 1', 'Principle Component 2']], ax, c, means_df, names[n])
			#plt.clf()
				all_pred.append(pred)
				print(acc)
				all_accs.append(acc)
			else:
				all_accs.append(0)
		else:
			all_accs.append(0)
	info['accuracy'] = all_accs
	final = pd.merge(info, driver_df['pathway'], left_index = True, right_index = True)
	final = final[(final['accuracy'] != 0)]
	final['classes'] = all_pred
	real_drivers = list(set(driver_names).intersection(preds.columns.values))
