#-*- coding: UTF-8 -*- 
import codecs
import os
import sys
import re

DATA_FORMAT = "(float)"


class Formula:
	def __init__(self, left, right):
		self.left = left
		self.right = right
	def all(self):
		return self.left + ' = ' + self.right


class Equation:
	def __init__(self):
		self.coefficient = []
		self.variable = []
		self.rightItem = ''

	def total(self, map):
		self.n = len(self.variable)
		back = "\n" + "	FEI newFEI(" + str(self.n)+", " + map + ");\n"
		for i in range(self.n):
			back += "	newFEI.set("+self.variable[i] + ", " +self.coefficient[i] + ");\n"
		back += "	newFEI.rightItem = " + self.rightItem + ";\n"
		return back


def findLeftCurlyBrace(str, p):
	# 从str的p位置开始向后寻找匹配的‘{’，并返回其位置
	curlyBraceCount = 0
	while str[p] != '{' or curlyBraceCount != 0:
		if str[p] == '{':
			curlyBraceCount += 1
		elif str[p] == '}':
			curlyBraceCount -= 1
		p -= 1
	return p


def findRightCurlyBrace(str, p):
	# 从str的p位置开始向后寻找匹配的‘}’，并返回其位置
	curlyBraceCount = 0
	while str[p] != '}' or curlyBraceCount != 0:
		if str[p] == '{':
			curlyBraceCount += 1
		elif str[p] == '}':
			curlyBraceCount -= 1
		p += 1
	return p

greekLetter = [ r'\\alpha',
				r'\\beta',
				r'\\gamma',
				r'\\delta',
				r'\\epsilon',
				r'\\zeta',
				r'\\eta',
				r'\\theta',
				r'\\iota',
				r'\\kappa',
				r'\\lambda',
				r'\\mu',
				r'\\xi',
				r'\\pi',
				r'\\rho',
				r'\\sigma',
				r'\\varsigma',
				r'\\upsilon',
				r'\\phi',
				r'\\chi',
				r'\\psi',
				r'\\omega',
				r'\\Gamma',
				r'\\Delta',
				r'\\varepsilon',
				r'\\Theta',
				r'\\varkappa',
				r'\\Lambda',
				r'\\nu',
				r'\\Xi',
				r'\\Pi',
				r'\\varrho',
				r'\\Sigma',
				r'\\tau',
				r'\\Upsilon',
				r'\\Phi',
				r'\\Psi',
				r'\\Omega',
				r'\\digamma',
				r'\\vartheta',
				r'\\varpi',
				r'\\varsigma',
				r'\\varphi']
def latexTranslate(str):
	# 指数处理
	i = re.search(r'\}\^(?P<num>\d+)', str)
	while i :
		p = findLeftCurlyBrace(str, i.span()[0]-1)
		str = str[:p] + 'pow(' + str[p+1:i.span()[0]] + ', ' + i.group('num') + ')' + str[i.span()[1]:]
		i = re.search(r'\}\^(?P<num>\d+)', str)
	# 转换下标与数组访问
	i = re.search(r'[A-Za-z_]\w*_\{', str)
	while i :
		start = i.span()[0]
		end = findRightCurlyBrace(str, i.span()[1])+1
		subStr = str[start:end]
		subStr = re.sub(r'_\{', '[', subStr)
		subStr = re.sub(r',', '][', subStr)
		subStr = re.sub(r'\}', ']', subStr)
		str = str[:start] + subStr + str[end:]
		i = re.search(r'[A-Za-z_]\w*_\{', str)
	# 转换分数
	i = re.search(r'\\frac\{', str)
	while i :
		p0 = i.span()[1]#位于第一个{后的字符
		p1 = findRightCurlyBrace(str, p0)#位于第一个}
		p2 = findRightCurlyBrace(str, p1+2)#位于第二个}
		str = str[:i.span()[0]] + '(' + str[p0:p1] + r')/(' + str[p1+2:p2] + ')' + str[p2+1:]
		i = re.search(r'\\frac\{', str)
	# 去除转义字符*
	str = re.sub(r'\\\*', '*', str)
	# 去除希腊字母和其他的转义字符
	for i in greekLetter:
		str = re.sub(i, i[2:], str)
		str = re.sub(i[:2]+'\{'+i[2:]+'\}', i[2:], str)
	# 添加忽略的乘法
	i = re.search(r'(^|[^A-Za-z])(?P<num>[\d|.]+)[A-Za-z\(]', str)
	while i :
		str = re.sub(i.group('num'), i.group('num')+'*', str)
		i = re.search(r'(^|[^A-Za-z])(?P<num>\d+)[A-Za-z\(]', str)
	# 处理整数精度转化
	i = re.search(r'(^|[^A-Za-z\)])(?P<num>\d+)(\n|$|[^A-Za-z])', str)
	while i :
		str = re.sub(i.group('num'), DATA_FORMAT + i.group('num'), str)
		i = re.search(r'(^|[^A-Za-z)])(?P<num>\d+)(\n|$|[^A-Za-z])', str)
	return str


class FEI:
	'''
	file : string[] FEI文件每行的内容
	formula : class Formula{} 算式数组
	equation : class Equation{} 方程数组
	'''
	def openFile(self):
		try:
			filePointer = codecs.open(sys.path[0] + r"/../doc/FEI.md","r","utf-8")
		except Exception:
			print("FEI.md doesn't exist in doc!")
			sys.exit(1)
		self.file = filePointer.readlines()
		filePointer.close()

	def recognizeFormula(self, name, statement):
		# 删除####
		name = re.sub(r'####', '', name)
		# 删除注释
		name = re.sub(r':.*', '', name)
		# 删除空格及回车
		name = re.sub(r'\s', '', name)
		# 删除$
		statement = re.sub(r'\$', '', statement)
		# 删除空格
		statement = re.sub(r'\s', '', statement)
		#print(name);#print(statement)

		# 分解等号左右两边并分别处理
		seperate = re.match(r'(?P<left>.+)=(?P<right>.+)', statement)
		left = seperate.group('left')
		right = seperate.group('right')
		left = latexTranslate(left)
		right = latexTranslate(right)
		#print(left, right)
		formula = Formula(left, right)
		self.formula[name] = formula
	
	def recognizeEquation(self, name, statement):
		#删除####@
		name = re.sub('####@', '', name)
		#删除注释
		name = re.sub(r':.*', '', name)
		#删除空格及回车
		name = re.sub(r'\s', '', name)
		#删除$
		statement = re.sub(r'\$', '', statement)
		#删除空格
		statement = re.sub(r'\s', '', statement)
		#print(name);print(statement)

		seperate = re.match(r'(?P<left>.+)=(?P<right>.+)', statement)
		left = seperate.group('left')
		right = seperate.group('right')
		#print(left, right)
		equation = Equation()
		equation.rightItem = latexTranslate(right)
		#print(equation.rightItem)
		start = 0
		for i in re.finditer(r'\\underline{', left):
			#print(i)
			if start == i.span()[0]:
				equation.coefficient.append("1")
			elif start == i.span()[0]-1:
				equation.coefficient.append(left[start:i.span()[0]] + "1")
			else:
				equation.coefficient.append(latexTranslate(left[start:i.span()[0]]))
			p = findRightCurlyBrace(left, i.span()[1])
			subStr = left[i.span()[1]:p]
			subStr = re.sub(r'_\{', 'Map(', subStr)
			subStr = re.sub(r'\}', ')', subStr)
			equation.variable.append(subStr)
			start = p+1
		#print(equation.variable)
		#print(equation.coefficient)
		self.equation[name] = equation

	def __init__(self):
		self.openFile()
		self.formula = {}
		self.equation = {}
		for i in range(len(self.file)):
			if self.file[i].startswith("####@"):
				if i+1 >= len(self.file):
					print("Rules not complete at the end of FEI.md !")
					sys.exit(1)
				try:
					self.recognizeEquation(self.file[i], self.file[i+1])
				except Exception:
					print('Error at line', i+1, 'or', i+2)
					sys.exit(1)
			elif self.file[i].startswith("####"):
				if i+1 >= len(self.file):
					print("Rules not complete at the end of FEI.md !")
					sys.exit(1)
				try:
					self.recognizeFormula(self.file[i], self.file[i+1])
				except Exception:
					print('Error at line', i+1, 'or', i+2)
					sys.exit(1)


def compileCpp(inputFileName, outputFileName):
	try:
		inputFile = codecs.open(inputFileName, "r", "utf-8")
	except Exception:
		print("Can not open" + inputFileName)
		sys.exit(1)
	try:
		outputFile = codecs.open(outputFileName, "w", "utf-8")
	except Exception:
		print("Can not create", outputFileName)
		sys.exit(1)
	inputFileData = inputFile.readlines()
	for i in inputFileData:
		target = re.search(r'##\s*(?P<name>\w+)', i)
		if target:
			targetName = target.group('name')
			if targetName in FEI.formula:
				if re.search(r'(=|return).*##\s*\w+', i):
					newString = re.sub(r'##.*\w+', FEI.formula[targetName].right, i)
					outputFile.write(newString)
				else:
					newString = re.sub(r'##.*\w+', FEI.formula[targetName].all(), i)
					outputFile.write(newString)
			elif target.group('name') in FEI.equation:
				map = re.search(r'\{(?P<map>.+)\}', i)
				map = map.group('map')
				newString = re.sub(r'##.*', FEI.equation[targetName].total(map), i)
				outputFile.write(newString)
			else:
				print("Can not find rules in FEI.md by name", targetName)
				sys.exit(1)
		else:
			outputFile.write(i)

	inputFile.close()
	outputFile.close()

if __name__ == "__main__":
	FEI = FEI()
	'''files = os.listdir()
	for i in files:
		if re.match(r'.*\.cpp', i):
			compileCpp(i)'''
	if len(sys.argv) < 2:
		print('No file is indicated!')
		sys.exit(1)
	inputFileName = sys.argv[1]
	outputFileName = sys.argv[2];
	compileCpp(inputFileName, outputFileName)
	print('FEI compile', inputFileName, 'finished!')

