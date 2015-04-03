from whatshap.priorityqueue import PriorityQueue

def test_queue():
	pq = PriorityQueue()
	pq.push(10, ord('a'))
	pq.push(5, ord('b'))
	pq.push(12, ord('c'))
	pq.push(3, ord('d'))
	assert len(pq) == 4
	assert pq.pop() == (12, ord('c'))
	assert pq.pop() == (10, ord('a'))
	assert pq.pop() == (5, ord('b'))
	assert pq.pop() == (3, ord('d'))

def test_queue2():
	pq = PriorityQueue()
	pq.push(1, ord('a'))
	pq.push(50, ord('b'))
	pq.push(2, ord('c'))
	pq.push(30, ord('d'))
	assert len(pq) == 4
	assert pq.pop() == (50, ord('b'))
	assert pq.pop() == (30, ord('d'))
	assert pq.pop() == (2, ord('c'))
	assert pq.pop() == (1, ord('a'))


def test_change_score():
	pq = PriorityQueue()
	pq.push(10, ord('a'))
	pq.push(5, ord('b'))
	pq.change_score(ord('a'), 2)
	pq.push(12, ord('c'))
	pq.push(3, ord('d'))
	pq.change_score(ord('c'), 1)
	pq.change_score(ord('d'), 15)
	assert len(pq) == 4
	assert pq.pop() == (15, ord('d'))
	assert pq.pop() == (5, ord('b'))
	assert pq.pop() == (2, ord('a'))
	assert pq.pop() == (1, ord('c'))

def test_change_score_sorting():
	pq=PriorityQueue()
	pq.push(50,ord('1'))
	pq.push(40,ord('2'))
	pq.push(30,ord('3'))
	pq.push(20,ord('4'))
	pq.push(10,ord('5'))
	pq.change_score(ord('5'),100)
	pq.change_score(ord('2'),45)
	assert len(pq) == 5
	assert pq.pop() == (100,ord('5'))
	#assert pq.pop() == (45,ord('2'))
	pq.push(60,ord('8'))
	assert pq.pop() == (60,ord('8'))
	pq.change_score(ord('2'),40)
	assert pq.pop()== (50,ord('1'))
	assert pq.pop() == (40,ord('2'))

def test_is_empty():
	pq=PriorityQueue()
	assert pq.is_empty()
	pq.push(10,ord('A'))
	assert not pq.is_empty()
	pq.pop()
	assert pq.is_empty()
	pq.push(9,ord('B'))
	assert not pq.is_empty()
	pq.push(8,ord('C'))
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert pq.is_empty()
	pq.push(7,ord('D'))
	assert not pq.is_empty()
	pq.push(6,ord('E'))
	assert not pq.is_empty()
	pq.push(5,ord('F'))
	assert not pq.is_empty()
	pq.push(4,ord('G'))
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert pq.is_empty()


def test_tuple_score():
	pq= PriorityQueue()
	assert pq.is_empty()
	pq.push((4,0,2),ord('A'))
	(score,item ) = pq.pop()
	assert score== (4,0,2)
	assert item==ord('A')
	pq.push((1,0,2),ord('B'))
	pq.push((2,0,2),ord('C'))
	pq.push((3,0,2),ord('D'))
	pq.push((40,0,2),ord('E'))
	pq.push((0,0,2),ord('F'))
	pq.push((50,0,2),ord('G'))
	(score,item )=pq.pop()
	assert score == (50,0,2)
	assert item ==ord('G')
	(score,item )=pq.pop()
	assert score == (40,0,2)
	assert item ==ord('E')
	(score,item )=pq.pop()
	assert score == (3,0,2)
	assert item ==ord('D')
	(score,item )=pq.pop()
	assert score == (2,0,2)
	assert item ==ord('C')
	pq.pop()
	pq.pop()
	assert pq.is_empty()

def test_tuple_score_sorting():
	pq= PriorityQueue()
	pq.push((10,0,0),ord('B'))
	pq.push((10,2,6),ord('C'))
	pq.push((10,3,2),ord('D'))
	pq.push((10,4,3),ord('E'))
	pq.push((10,2,2),ord('F'))
	pq.push((10,0,2),ord('G'))
	(score,item ) = pq.pop()
	assert score== (10,4,3)
	assert item==ord('E')
	(score,item ) = pq.pop()
	assert score== (10,3,2)
	assert item==ord('D')
	(score,item ) = pq.pop()
	assert score== (10,2,6)
	assert item==ord('C')
	(score,item ) = pq.pop()
	assert score== (10,2,2)
	assert item==ord('F')
	(score,item ) = pq.pop()
	assert score== (10,0,2)
	assert item==ord('G')
	pq.push((1,10,4),ord('X'))
	pq.push((5,0,6),ord('Y'))
	pq.push((1,8,2),ord('Z'))
	pq.change_score(ord('Y'),(100,100,100))
	pq.change_score(ord('Z'),(0,0,0))
	(score,item ) = pq.pop()
	assert score== (100,100,100)
	assert item==ord('Y')
	(score,item ) = pq.pop()
	assert score== (10,0,0)
	assert item==ord('B')
	(score,item ) = pq.pop()
	assert score== (1,10,4)
	assert item==ord('X')
	(score,item ) = pq.pop()
	assert score== (0,0,0)
	assert item==ord('Z')
	assert pq.is_empty()
