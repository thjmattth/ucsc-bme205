	def traceit(algorithm):
		"""
		performs traceback on local alignment

		:return str, master-slave alignment in A2M
		"""

		#row, col is the current position of the traceback in the matrices 
		#self.M, self.IR, and self.IC.
		
		row = len(self.row_seq)
		col = len(self.col_seq)

		# records traceback as a list aligned bases in reverse
		trace_list = []

		if algorithm = 'local':
			# appends insertions/deletions to best alignment in the local case
			while(row > self.best_row):
				trace_list += '-'
				row -= 1
			while(col > self.best_col):
				trace_list += self.col_seq[col-1].lower()
				col -= 1

			hi_score = self.M[self.best_row][self.best_col]
			hi_score_neighbor = 0
			end_align = False
			
		
		if algorithm = 'global':
			# appends insertions/deletions to best alignment in the local case
			hi_score_neighbor, hi_score = max(enumerate([self.M[row][col], self.IR[row][col], self.IC[row][col]]), key=lambda x:x[1])
			
			
		while((row != 0 or col != 0) and not end_align):

            if hi_score_neighbor == 0:
                trace_list += self.col_seq[col-1].upper()
                subst = self.subst[(self.row_seq[row-1], self.col_seq[col-1])]

                row -= 1
                col -= 1

                if(hi_score == subst):
                    end_align = True
                elif(hi_score == subst + self.M[row][col]):
                    hi_score_neighbor = 0
                    hi_score = self.M[row][col]
                elif(hi_score == subst + self.IR[row][col]):
                    hi_score_neighbor = 1
                    hi_score = self.IR[row][col]
                elif(hi_score == subst + self.IC[row][col]):
                    hi_score_neighbor = 2
                    hi_score = self.IC[row][col]

            elif hi_score_neighbor == 1:
                trace_list += '-'

                row -= 1

                if(hi_score + self.start == self.M[row][col]):
                    hi_score_neighbor = 0
                    hi_score = self.M[row][col]
                elif(hi_score + self.extend == self.IR[row][col]):
                    hi_score_neighbor = 1
                    hi_score = self.IR[row][col]
                elif(hi_score + self.double == self.IC[row][col]):
                    hi_score_neighbor = 2
                    hi_score = self.IC[row][col]

            elif hi_score_neighbor == 2:
                trace_list += self.col_seq[col-1].lower()

                col -= 1

                if(hi_score + self.start == self.M[row][col]):
                    hi_score_neighbor = 0
                    hi_score = self.M[row][col]
                elif(hi_score + self.double == self.IR[row][col]):
                    hi_score_neighbor = 1
                    hi_score = self.IR[row][col]
                elif(hi_score + self.extend == self.IC[row][col]):
                    hi_score_neighbor = 2
                    hi_score = self.IC[row][col]

        while(col > 0):
            trace_list += self.col_seq[col-1].lower()
            col -= 1
        while(row > 0):
            trace_list += '-'
            row -= 1

        print('this is rev_alignment list')
        print(trace_list)

        return ''.join(trace_list[::-1])
