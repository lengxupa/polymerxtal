class Person:
    def __init__(self, name, surname):
        self.name = name
        self.surname = surname
        self.id = self.generate_id()
    
    def generate_id(self):
        id_hash = 0
        for s in self.name:
            id_hash += ord(s)
        for s in self.surname:
            id_hash *= ord(s)
        return id_hash % 1000000000
    
    def __str__(self):
        return f'{self.surname}, {self.name}\tID: {self.id}'


class Student(Person):
    def __init__(self, name, surname):
        self.courses = []
        super().__init__(name, surname)
    
    def __str__(self):
        return super().__str__() + f'\nCourses:\n{self.courses}'
        
    def enroll(self, new_course):
        self.courses.append(new_course)
        
    def drop_course(self, course):
        self.courses.remove(course)


class Faculty(Person):
    def __init__(self, name, surname, position, salary):
        self.position = position
        self.salary = salary
        self.courses = []
        super().__init__(name, surname)
 
    def __str__(self):
        return super().__str__() + f'\nCourses:\n{self.courses}'
 
    def assign_course(self, new_course):
        self.courses.append(new_course)
 
    def unassign_course(self, course):
        self.courses.remove(course)

